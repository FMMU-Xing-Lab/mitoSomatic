#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

# 从我们刚刚写的两个脚本中引入核心类
from mut_features_extractor import MutationFeatureExtractor
from rf_classifier import MutationPredictor

def run_pipeline(input_file, output_dir):
    # 1. 检查输入文件与创建输出目录
    if not os.path.exists(input_file):
        print(f"【致命错误】找不到输入文件: {input_file}")
        sys.exit(1)
        
    if not os.path.exists(output_dir):
        print(f"检测到输出目录不存在，正在创建: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)

    # 2. 定义中间特征文件和最终预测结果文件的路径
    base_name = os.path.basename(input_file)
    features_outfile = os.path.join(output_dir, f"{base_name}.features.txt")
    final_outfile = os.path.join(output_dir, f"{base_name}.predicted.txt")

    print("\n" + "="*50)
    print("      启动线粒体突变分析管线 (Pipeline)")
    print("="*50)

    # ==========================
    # 步骤一：提取9个多维特征
    # ==========================
    print("\n[步骤 1/2] 开始提取多维特征...")
    extractor = MutationFeatureExtractor()
    extractor.process_file(input_file, features_outfile)

    # ==========================
    # 步骤二：利用随机森林进行分类预测
    # ==========================
    print("\n[步骤 2/2] 开始加载随机森林模型并预测...")
    predictor = MutationPredictor()
    predictor.predict(features_outfile, final_outfile)

    print("\n" + "="*50)
    print(f"🎉 全部流程执行成功！")
    print(f"中间特征文件保存于: {features_outfile}")
    print(f"最终预测结果保存于: {final_outfile}")
    print("="*50 + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("【错误】参数不正确！")
        print("【用法】: python main_pipeline.py <输入突变文件> <输出结果路径>")
        print("【示例】: python main_pipeline.py my_input.txt ./result_folder/")
        sys.exit(1)
        
    input_txt = sys.argv[1]
    out_dir = sys.argv[2]
    
    run_pipeline(input_txt, out_dir)
