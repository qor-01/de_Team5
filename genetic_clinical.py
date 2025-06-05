import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.font_manager as fm
import os

# 한글 폰트 설정 (말굽체)
font_path = 'C:/Windows/Fonts/malgun.ttf'
plt.rcParams['font.family'] = fm.FontProperties(fname=font_path).get_name()

# chrY 파일 불러오기
df = pd.read_csv('KOVA_chrY_p.csv')

# ClinVar 복수 표기에서 첫 항목만 추출
df['clinvar_CLNSIG_simple'] = df['clinvar_CLNSIG'].str.split(',').str[0].str.strip()

# 병원성 여부 판별 함수
def is_pathogenic(clnsig):
    return any(x in clnsig for x in ['Pathogenic', 'Likely_pathogenic'])

df['clinvar_pathogenic'] = df['clinvar_CLNSIG_simple'].apply(is_pathogenic)

# 숫자로 변환
df['gwas_pvalue'] = pd.to_numeric(df['gwas_pvalue'], errors='coerce')

# 유의미한 GWAS 변이 추출
gwas_significant = df[df['gwas_pvalue'] < 5e-8]

# 분포 시각화
plt.figure(figsize=(10, 5))
sns.histplot(df['gwas_pvalue'], bins=50, log_scale=True, color='purple')
plt.axvline(5e-8, color='red', linestyle='--', label='Genome-wide significance')
plt.title('GWAS p-value 분포')
plt.xlabel('p-value (log scale)')
plt.ylabel('Count')
plt.legend()
plt.tight_layout()
plt.show()

# GWAS 위험요소 클래스 시각화
plt.figure(figsize=(8, 4))
sns.countplot(data=gwas_significant, y='gwas_class', palette='coolwarm')
plt.title('GWAS 연관 클래스 (유의한 변이)')
plt.xlabel('Count')
plt.ylabel('GWAS Class')
plt.tight_layout()
plt.show()

# 숫자 변환
df['gtex_qval'] = pd.to_numeric(df['gtex_qval'], errors='coerce')

# 유의미한 GTEx 변이 추출
gtex_significant = df[df['gtex_qval'] < 0.05]

# 시각화
plt.figure(figsize=(8, 4))
sns.countplot(data=gtex_significant, y='gtex_class', palette='Set2')
plt.title('GTEx 유전자 발현 영향 (q < 0.05)')
plt.xlabel('Count')
plt.ylabel('GTEx Class')
plt.tight_layout()
plt.show()

combined_filter = (
    df['clinvar_pathogenic'] &
    (df['gwas_pvalue'] < 5e-8) &
    (df['gtex_qval'] < 0.05)
)

important_variants = df[combined_filter]

# 결과 확인
print(f"임상적으로 중요한 변이 수: {important_variants.shape[0]}")
print(important_variants[['clinvar_CLNSIG', 'gwas_class', 'gtex_class', 'gwas_pvalue', 'gtex_qval']].head())

summary = []
chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

for chrom in chromosomes:
    try:
        file_path = f'KOVA_{chrom}_p.csv'
        if not os.path.exists(file_path):
            print(f"📁 파일 없음: {file_path}")
            summary.append({'Chromosome': chrom, 'ImportantVariants': 0})
            continue

        df = pd.read_csv(file_path)

        # ClinVar 병원성 처리
        df['clinvar_CLNSIG_simple'] = df['clinvar_CLNSIG'].str.split(',').str[0].str.strip()
        df['clinvar_pathogenic'] = df['clinvar_CLNSIG_simple'].apply(is_pathogenic)

        # GWAS, GTEx 숫자 변환
        df['gwas_pvalue'] = pd.to_numeric(df['gwas_pvalue'], errors='coerce')
        df['gtex_qval'] = pd.to_numeric(df['gtex_qval'], errors='coerce')

        # 필터링
        filter_all = (
            df['clinvar_pathogenic'] &
            (df['gwas_pvalue'] < 5e-8) &
            (df['gtex_qval'] < 0.05)
        )

        count = df[filter_all].shape[0]
        summary.append({'Chromosome': chrom, 'ImportantVariants': count})
        print(f"{chrom}: {count}개의 중요한 변이 발견")

    except Exception as e:
        print(f"❌ {chrom} 처리 중 오류 발생: {e}")
        summary.append({'Chromosome': chrom, 'ImportantVariants': 0})

# 결과 요약
summary_df = pd.DataFrame(summary)
print("\n최종 요약:")
print(summary_df)

if summary_df['ImportantVariants'].sum() == 0:
    print("⚠️ 모든 염색체에서 중요한 변이가 발견되지 않았습니다. 그래프를 그릴 수 없습니다.")
else:
    plt.figure(figsize=(12, 6))
    sns.barplot(data=summary_df, x='Chromosome', y='ImportantVariants', palette='viridis')
    plt.title('염색체별 임상적으로 중요한 변이 수')
    plt.xlabel('Chromosome')
    plt.ylabel('중요 변이 수')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
