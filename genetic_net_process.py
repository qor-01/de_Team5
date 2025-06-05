import pandas as pd
import numpy as np

df = pd.read_csv("KOVA_chr1_p.csv", low_memory=False)

# kova_AF >= 0.01 데이터만
df_af_filtered = df[df['kova_AF'] >= 0.01]
df_af_filtered.head(10)

df_af_filtered.shape

# 예측 정보 없는 행 삭제
mask = (
    (df_af_filtered['clinvar_pathogenic'] == 0) &
    (df_af_filtered['gwas_class'] == 0) &
    (df_af_filtered[['CADD', 'ReMM', 'FunSeq2', 'LINSIGHT']].isnull().all(axis=1))
)

# 제거
df_filtered = df_af_filtered[~mask]

print(df)
print(df_filtered.shape)

# 파일 저장
df_filtered.to_csv('KOVA_chrY_net_p.csv', index=False, encoding='utf-8')