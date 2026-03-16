#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess

# ================= 绝对路径配置区 =================
# 建议这里直接写死您的绝对路径
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# ANNOVAR 脚本与数据库配置
ANNOTATE_VAR_PL = os.path.join(BASE_DIR, "annotate_variation.pl")
HUMANDB_DIR = os.path.join(BASE_DIR, "humandb")
BUILD_VER = "hg19"
DB_TYPE = "MT_ensGene"
DEFAULT_CHR = "chrM"

# 8个自定义特征来源文件配置
REFBASE_FILE = os.path.join(BASE_DIR, "refbase.txt")
MUT_DIST_FILE = os.path.join(BASE_DIR, "Mutational_distribution.txt")
MUT_ASSESSOR_FILE = os.path.join(BASE_DIR, "Mutation assessor.txt")
DBSNP_FILE = os.path.join(BASE_DIR, "dbSNP-fix.txt")
VAF_MITOMAP_FILE = os.path.join(BASE_DIR, "VAF in mitomap.txt")
PHYLOTREE_FILE = os.path.join(BASE_DIR, "Phylotree haplogroup status.txt")
MTDB_FILE = os.path.join(BASE_DIR, "mtDB polymorphic sites.txt")
# =================================================


class MutationFeatureExtractor:
    """
    线粒体突变多维特征提取工具包
    """
    def __init__(self):
        self.check_files()
        print("正在加载本地数据库到内存，请稍候...")
        self.load_databases()

    def check_files(self):
        """检查所需文件是否存在"""
        files_to_check = [
            ANNOTATE_VAR_PL, REFBASE_FILE, MUT_DIST_FILE, 
            MUT_ASSESSOR_FILE, DBSNP_FILE, VAF_MITOMAP_FILE, 
            PHYLOTREE_FILE, MTDB_FILE
        ]
        for f in files_to_check:
            if not os.path.exists(f):
                raise FileNotFoundError(f"找不到必要文件: {f}")

    def load_databases(self):
        """加载所有特征映射表到内存 (使用字典/集合提升查询速度)"""
        
        # 1. refbase (位点 -> 64种背景 Context)
        self.refbase_dict = {}
        with open(REFBASE_FILE, 'r') as f:
            for line in f:
                if line.startswith("site"): continue
                parts = line.strip().split()
                if len(parts) >= 3:
                    self.refbase_dict[parts[0]] = parts[2]

        # 2. Mutational_distribution (位点 -> PRcoding/Noncoding 等)
        self.mut_dist_dict = {}
        with open(MUT_DIST_FILE, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4:
                    self.mut_dist_dict[parts[0]] = parts[3]

        # 3. Mutation_assessor (位点Ref>Alt -> 影响程度)
        self.mut_assessor_dict = {}
        with open(MUT_ASSESSOR_FILE, 'r') as f:
            for line in f:
                parts = line.strip().split('\t') # 使用tab分割以防最后列为空
                if len(parts) >= 1:
                    mut_key = parts[0]
                    # 如果列数足够且第五列有内容，则写入，否则记为NAN
                    if len(parts) >= 5 and parts[4].strip() != "":
                        self.mut_assessor_dict[mut_key] = parts[4].strip()
                    else:
                        self.mut_assessor_dict[mut_key] = "NAN"

        # 4. dbSNP (集合: 位点Ref>Alt)
        self.dbsnp_set = set()
        with open(DBSNP_FILE, 'r') as f:
            for line in f:
                if line.startswith("index"): continue
                self.dbsnp_set.add(line.strip())

        # 5. VAF in mitomap (位点Ref>Alt -> VAF值)
        self.vaf_mitomap_dict = {}
        with open(VAF_MITOMAP_FILE, 'r') as f:
            for line in f:
                if line.startswith("index"): continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    self.vaf_mitomap_dict[parts[0]] = parts[1]

        # 6. Phylotree (集合: 位点Ref>Alt)
        self.phylotree_set = set()
        with open(PHYLOTREE_FILE, 'r') as f:
            for line in f:
                if line.startswith("pos"): continue
                self.phylotree_set.add(line.strip())

        # 7. mtDB (集合: 位点Ref>Alt)
        self.mtdb_set = set()
        with open(MTDB_FILE, 'r') as f:
            for line in f:
                if line.startswith("pos"): continue
                self.mtdb_set.add(line.strip())

        print("本地数据库加载完成！")

    def _map_annovar_func(self, annovar_func):
        func_lower = annovar_func.lower()
        if 'nonsynonymous' in func_lower: return 'Nonsynonymous'
        elif 'synonymous' in func_lower: return 'synonymous'
        elif 'stopgain' in func_lower: return 'stopgain'
        elif 'stoploss' in func_lower: return 'stoploss'
        else: return 'Nan'

    def run_annovar(self, input_file):
        """运行ANNOVAR并返回字典映射: key='pos_ref_alt', value='Annotation'"""
        temp_avinput = input_file + ".temp.avinput"
        out_prefix = input_file + ".annovar_out"

        print("正在调用 ANNOVAR 获取注释特征...")
        with open(input_file, 'r') as fin, open(temp_avinput, 'w') as fout:
            for line in fin:
                parts = line.strip().split()
                if len(parts) >= 4:  # Sample, Pos, Ref, Alt, VAF
                    pos, ref, alt = parts[1], parts[2], parts[3]
                    fout.write(f"{DEFAULT_CHR}\t{pos}\t{pos}\t{ref}\t{alt}\n")

        cmd = [
            "perl", ANNOTATE_VAR_PL, "-geneanno",
            "-buildver", BUILD_VER, "-dbtype", DB_TYPE,
            temp_avinput, HUMANDB_DIR, "-out", out_prefix
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

        anno_dict = {}
        exonic_file = f"{out_prefix}.exonic_variant_function"
        if os.path.exists(exonic_file):
            with open(exonic_file, 'r') as f_exonic:
                for line in f_exonic:
                    parts = line.strip().split('\t')
                    if len(parts) >= 8:
                        mut_type, start, ref, alt = parts[1], parts[4], parts[6], parts[7]
                        anno_dict[f"{start}_{ref}_{alt}"] = self._map_annovar_func(mut_type)

        # 清理临时文件
        for ext in [".temp.avinput", ".annovar_out.variant_function", ".annovar_out.exonic_variant_function", ".annovar_out.log"]:
            tmp_f = input_file + ext
            if os.path.exists(tmp_f): os.remove(tmp_f)

        return anno_dict

    def process_file(self, input_file, output_file):
        """处理输入文件，输出包含所有特征的表格"""
        
        # 先提取 ANNOVAR 字典
        annovar_dict = self.run_annovar(input_file)

        print("正在进行多维特征匹配与融合...")
        with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
            # 写入表头（包含您需要的9个新特征）
            headers = [
                "Sample", "Position", "Ref", "Alt", "VAF",
                "ANNOVAR", "sub_type", "refbase", "Mutational_distribution",
                "Mutation_assessor", "dbSNP", "VAF_mitomap", "Phylotree", "mtDB"
            ]
            fout.write("\t".join(headers) + "\n")

            for line in fin:
                line = line.strip()
                if not line or line.startswith("#"): continue
                parts = line.split()
                if len(parts) < 5:
                    print(f"警告：跳过列数不足的行 -> {line}")
                    continue

                sample = parts[0]
                pos = parts[1]
                ref = parts[2]
                alt = parts[3]
                vaf = parts[4]

                # 构建联合主键
                annovar_key = f"{pos}_{ref}_{alt}"
                mut_str_key = f"{pos}{ref}>{alt}"  # e.g. "3469C>T"
                sub_type = f"{ref}>{alt}"          # e.g. "C>T"

                # === 特征提取逻辑 ===
                # 1. ANNOVAR
                feat_annovar = annovar_dict.get(annovar_key, "Nan")
                # 2. sub_type
                feat_sub_type = sub_type
                # 3. refbase (查询第三列)
                feat_refbase = self.refbase_dict.get(pos, "NAN")
                # 4. Mutational_distribution (查询第四列)
                feat_mut_dist = self.mut_dist_dict.get(pos, "NAN")
                # 5. Mutation_assessor
                feat_mut_assessor = self.mut_assessor_dict.get(mut_str_key, "NAN")
                # 6. dbSNP
                feat_dbsnp = "1" if mut_str_key in self.dbsnp_set else "0"
                # 7. VAF_mitomap
                feat_vaf_mitomap = self.vaf_mitomap_dict.get(mut_str_key, "NAN")
                # 8. Phylotree
                feat_phylotree = "1" if mut_str_key in self.phylotree_set else "0"
                # 9. mtDB
                feat_mtdb = "1" if mut_str_key in self.mtdb_set else "0"

                # 写入汇总结果
                out_row = [
                    sample, pos, ref, alt, vaf,
                    feat_annovar, feat_sub_type, feat_refbase, feat_mut_dist,
                    feat_mut_assessor, feat_dbsnp, feat_vaf_mitomap, feat_phylotree, feat_mtdb
                ]
                fout.write("\t".join(out_row) + "\n")

        print(f"处理完毕！特征文件已成功保存至: {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("【用法】: python mut_features_extractor.py <输入突变文件> <输出特征文件>")
        print("【输入格式】: Sample  Position  Ref  Alt  VAF (用Tab或空格隔开，无表头)")
        sys.exit(1)
        
    in_txt = sys.argv[1]
    out_txt = sys.argv[2]
    
    # 实例化并运行流程
    extractor = MutationFeatureExtractor()
    extractor.process_file(in_txt, out_txt)
