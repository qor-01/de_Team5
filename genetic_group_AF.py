import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import glob
import os

# UMAP 설치 여부 확인
try:
    import umap.umap_ as umap
    use_umap = True
except ImportError:
    print("UMAP이 없어 PCA만 사용합니다.")
    use_umap = False

# 한글 폰트 (윈도우 기준)
plt.rcParams['font.family'] = 'Malgun Gothic'
plt.rcParams['axes.unicode_minus'] = False

# 🔁 여러 CSV 병합 (메모리 절약형)
file_paths = glob.glob("processed_data/KOVA_chr*_p.csv")
df_list = []
for file in file_paths:
    df = pd.read_csv(file, low_memory=False)
    df_list.append(df[[
        'kova_AF', 'jpn_AF', 'gnomad_eas_AF', 'gnomad_eur_AF',
        'gnomad_all_AF', 'chimpanzee_AF'
    ]])

df = pd.concat(df_list, ignore_index=True)
df = df.apply(pd.to_numeric, errors='coerce')  # 숫자 아닌 값 처리

# ---------------- 1. AF 분포 겹침 ----------------
plt.figure(figsize=(10, 6))
for col in df.columns:
    sns.kdeplot(df[col].dropna(), label=col, fill=True, alpha=0.3)
plt.title('AF 분포 비교 (KDE)')
plt.xlabel('Allele Frequency')
plt.legend()
plt.tight_layout()
plt.show()

# ---------------- 2. AF 차이 ----------------
df['KOVA_JPN_diff'] = abs(df['kova_AF'] - df['jpn_AF'])
df['KOVA_EAS_diff'] = abs(df['kova_AF'] - df['gnomad_eas_AF'])
df['KOVA_EUR_diff'] = abs(df['kova_AF'] - df['gnomad_eur_AF'])
df['KOVA_ALL_diff'] = abs(df['kova_AF'] - df['gnomad_all_AF'])
df['KOVA_CHIMP_diff'] = abs(df['kova_AF'] - df['chimpanzee_AF'])

plt.figure(figsize=(10, 6))
sns.boxplot(data=df[[
    'KOVA_JPN_diff', 'KOVA_EAS_diff', 'KOVA_EUR_diff', 'KOVA_ALL_diff', 'KOVA_CHIMP_diff'
]])
plt.title('AF 차이 분포 (|KOVA - 비교집단|)')
plt.ylabel('AF Difference')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# ---------------- 3. 산점도 (y=x 일치도) ----------------
plt.figure(figsize=(16, 4))
for i, col in enumerate(['jpn_AF', 'gnomad_eas_AF', 'gnomad_eur_AF', 'gnomad_all_AF']):
    plt.subplot(1, 4, i + 1)
    sns.scatterplot(x=df['kova_AF'], y=df[col], alpha=0.2, s=10)
    plt.plot([0, 1], [0, 1], 'r--')  # y=x
    plt.title(f'KOVA vs {col}')
    plt.xlabel('kova_AF')
    plt.ylabel(col)
plt.tight_layout()
plt.show()

# ---------------- 4. PCA 및 UMAP ----------------
af_df = df[['kova_AF', 'jpn_AF', 'gnomad_eas_AF', 'gnomad_eur_AF', 'gnomad_all_AF', 'chimpanzee_AF']].dropna()
scaled = StandardScaler().fit_transform(af_df)

# PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled)
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])

plt.figure(figsize=(6, 5))
sns.scatterplot(x='PC1', y='PC2', data=pca_df, alpha=0.3)
plt.title('PCA로 본 집단 유사성')
plt.tight_layout()
plt.show()

# UMAP (선택사항)
if use_umap:
    reducer = umap.UMAP(n_components=2, random_state=42)
    umap_result = reducer.fit_transform(scaled)
    umap_df = pd.DataFrame(umap_result, columns=['UMAP1', 'UMAP2'])

    plt.figure(figsize=(6, 5))
    sns.scatterplot(x='UMAP1', y='UMAP2', data=umap_df, alpha=0.3)
    plt.title('UMAP으로 본 집단 유사성')
    plt.tight_layout()
    plt.show()
