import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 한글 폰트 설정
plt.rcParams['font.family'] = 'Malgun Gothic'
plt.rcParams['axes.unicode_minus'] = False

# 사용할 AC 컬럼 정의
ac_cols = ['kova_AC', 'jpn_AC', 'gnomad_eas_AC', 'gnomad_eur_AC', 'gnomad_all_AC', 'chimpanzee_AC']

# 대용량 CSV에서 필요한 열만 샘플링하여 읽기
df = pd.read_csv(
    "C:/Users/user/PycharmProjects/데이터엔지니어/processed_data/merged_KOVA.csv",
    usecols=ac_cols,
    engine='python',
    nrows=100000  # 최대 10만 행만 읽기 (메모리 절약)
)

# 숫자형 변환
df[ac_cols] = df[ac_cols].apply(pd.to_numeric, errors='coerce')

# KDE 분포 시각화
plt.figure(figsize=(10, 6))
for col in ac_cols:
    sns.kdeplot(df[col].dropna().sample(n=5000, random_state=42), label=col, fill=True, alpha=0.3)
plt.title('AC 분포 비교 (KDE)')
plt.xlabel('Allele Count')
plt.legend()
plt.tight_layout()
plt.show()

# 집단 간 AC 차이 계산
df['KOVA_JPN_AC_diff'] = abs(df['kova_AC'] - df['jpn_AC'])
df['KOVA_EAS_AC_diff'] = abs(df['kova_AC'] - df['gnomad_eas_AC'])
df['KOVA_EUR_AC_diff'] = abs(df['kova_AC'] - df['gnomad_eur_AC'])
df['KOVA_ALL_AC_diff'] = abs(df['kova_AC'] - df['gnomad_all_AC'])
df['KOVA_CHIMP_AC_diff'] = abs(df['kova_AC'] - df['chimpanzee_AC'])

ac_diff_cols = [col for col in df.columns if "_AC_diff" in col]

# 박스플롯
sampled_df = df[ac_diff_cols].dropna().sample(n=5000, random_state=42)
plt.figure(figsize=(10, 6))
sns.boxplot(data=sampled_df)
plt.title('AC 차이 분포 (|KOVA - 비교집단|)')
plt.ylabel('AC Difference')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# log 스케일 산점도
plt.figure(figsize=(15, 5))
for i, col in enumerate(ac_cols[1:5]):
    plt.subplot(1, 4, i + 1)
    sns.scatterplot(
        x=np.log2(df['kova_AC'].dropna().sample(n=5000, random_state=42) + 1),
        y=np.log2(df[col].dropna().sample(n=5000, random_state=42) + 1),
        alpha=0.3, s=10
    )
    plt.plot([0, 20], [0, 20], 'r--')
    plt.title(f'log2(KOVA_AC+1) vs log2({col}+1)')
    plt.xlabel('KOVA')
    plt.ylabel(col)
plt.tight_layout()
plt.show()

# MA Plot 함수 정의
def plot_ma(ac1, ac2, label1, label2):
    log_ac1 = np.log2(ac1 + 1)
    log_ac2 = np.log2(ac2 + 1)
    A = 0.5 * (log_ac1 + log_ac2)
    M = log_ac1 - log_ac2

    plt.figure(figsize=(6, 5))
    sns.scatterplot(x=A, y=M, alpha=0.3, s=10)
    plt.axhline(0, color='r', linestyle='--')
    plt.title(f'MA Plot: {label1} vs {label2}')
    plt.xlabel('A = Mean log2(AC)')
    plt.ylabel('M = log2(AC1) - log2(AC2)')
    plt.tight_layout()
    plt.show()

# MA plot 실행
sample_size = 5000
plot_ma(df['kova_AC'].dropna().sample(n=sample_size, random_state=42),
        df['jpn_AC'].dropna().sample(n=sample_size, random_state=42), 'KOVA', 'JPN')
plot_ma(df['kova_AC'].dropna().sample(n=sample_size, random_state=42),
        df['gnomad_eas_AC'].dropna().sample(n=sample_size, random_state=42), 'KOVA', 'gnomAD_EAS')
plot_ma(df['kova_AC'].dropna().sample(n=sample_size, random_state=42),
        df['gnomad_eur_AC'].dropna().sample(n=sample_size, random_state=42), 'KOVA', 'gnomAD_EUR')
plot_ma(df['kova_AC'].dropna().sample(n=sample_size, random_state=42),
        df['gnomad_all_AC'].dropna().sample(n=sample_size, random_state=42), 'KOVA', 'gnomAD_ALL')
plot_ma(df['kova_AC'].dropna().sample(n=sample_size, random_state=42),
        df['chimpanzee_AC'].dropna().sample(n=sample_size, random_state=42), 'KOVA', 'Chimpanzee')