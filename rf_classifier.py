#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

# 动态获取当前脚本所在目录
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TRAIN_FILE = os.path.join(BASE_DIR, "train.txt")

class MutationPredictor:
    def __init__(self):
        if not os.path.exists(TRAIN_FILE):
            raise FileNotFoundError(f"【错误】找不到训练集文件: {TRAIN_FILE}")
        
        self.model = None
        self._build_and_train_model()

    def _preprocess_data(self, df):
        """对特征进行格式清洗"""
        df_clean = df.copy()
        # VAF_mitomap 中的 'NAN' 替换为 0.0，并转为浮点数
        df_clean['VAF_mitomap'] = df_clean['VAF_mitomap'].replace('NAN', 0.0).astype(float)
        # 将二进制特征强转为整数
        for col in ['dbSNP', 'Phylotree', 'mtDB']:
            df_clean[col] = df_clean[col].astype(int)
        # VAF转浮点数
        df_clean['VAF'] = df_clean['VAF'].astype(float)
        return df_clean

    def _build_and_train_model(self):
        """构建管道并用 train.txt 训练模型"""
        print("正在加载训练集并训练随机森林模型...")
        train_df = pd.read_csv(TRAIN_FILE, sep='\t')
        train_df = self._preprocess_data(train_df)

        # 定义特征列
        self.numeric_features = ['VAF', 'dbSNP', 'VAF_mitomap', 'Phylotree', 'mtDB']
        self.categorical_features = ['ANNOVAR', 'sub_type', 'refbase', 'Mutational_distribution', 'Mutation_assessor']

        X_train = train_df[self.numeric_features + self.categorical_features]
        y_train = train_df['label']

        # 设置预处理管道
        preprocessor = ColumnTransformer(
            transformers=[
                ('num', SimpleImputer(strategy='constant', fill_value=0), self.numeric_features),
                # handle_unknown='ignore' 可以有效防止测试集中出现训练集没见过的 categorical 特征导致的报错
                ('cat', OneHotEncoder(handle_unknown='ignore', sparse_output=False), self.categorical_features)
            ])

        # 设置完整的随机森林预测管道
        self.model = Pipeline(steps=[
            ('preprocessor', preprocessor),
            ('classifier', RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1))
        ])

        # 训练模型
        self.model.fit(X_train, y_train)
        print("模型训练完成！")

    def predict(self, features_file, output_file):
        """对提取好特征的文件进行预测并输出"""
        print("正在对目标数据进行预测...")
        test_df = pd.read_csv(features_file, sep='\t')
        clean_test_df = self._preprocess_data(test_df)

        # 提取相应的特征进行预测
        X_test = clean_test_df[self.numeric_features + self.categorical_features]
        predictions = self.model.predict(X_test)

        # 将预测结果写入原始dataframe中
        test_df['label'] = predictions

        # 导出结果
        test_df.to_csv(output_file, sep='\t', index=False)
        print(f"预测完毕！带标签的最终结果已保存至: {output_file}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("【用法】: python rf_classifier.py <待预测特征文件> <输出文件>")
        sys.exit(1)
    
    predictor = MutationPredictor()
    predictor.predict(sys.argv[1], sys.argv[2])
