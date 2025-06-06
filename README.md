# de_Team5

# 프로젝트명 : 유전자 변이 네트워크 분석


# 팀원
- 백채원 (서울여대 데이터사이언스학과)
- 이재진 (서울여대 데이터사이언스학과)
- 최은혜 (서울여대 데이터사이언스학과)


# 프로젝트 설명
이 프로젝트는 유전체 변이 데이터를 수집하고, 집단 간 차이 및 임상적 유의성을 분석하며,  
네트워크 모델링을 통해 변이 간 상호작용을 시각화하는 것을 목표로 합니다.

또한, 기능 예측 점수를 바탕으로 핵심 변이를 식별함으로써  
유전체 변이 분석의 정밀도를 높이는 데 기여합니다.

---

# 프로젝트 과정
1. **도메인 학습**  
   - 주요 컬럼 선정 및 유전체 변이 이해  

2. **EDA (탐색적 데이터 분석)**  

3. **전처리**  
   - 결측치 처리
   - 네트워크 구현 (노드, 엣지)  

4. **임상적 유의성 분석**  

5. **집단 간 차이 분석**
   - AC
   - AF

6. **네트워크 모델링**  
   - 허브 노드 & 중심 노드 탐색  
   - 커뮤니티 탐지  
   - 중요도 측정  

---

# 데이터 출처
본 프로젝트는 **KOVA (Korean Variant Archive)**에서 제공하는 유전체 변이 데이터를 활용하였습니다.  
KOVA는 한국인을 대상으로 한 대규모 유전체 변이 데이터베이스로,  
집단 특이적 변이 분석 및 유전체 연구에 활용됩니다.  

[KOVA 공식 웹사이트](https://www.kobic.re.kr/kova/)에서 자세한 정보를 확인할 수 있습니다.

---

# 사용 기술 및 도구
- **언어**: Python  
- **라이브러리**: pandas, numpy, matplotlib, seaborn, scikit-learn (sklearn), glob, tqdm  
- **데이터베이스**: Neo4j  

---

# 주요 결과
- **희소한 네트워크 구조**: 유전자 간 직접 연결은 드물고, 전체적인 상호작용 밀도는 낮음  
- **기능적 모듈 존재**: 일부 유전자들은 작고 독립적인 클러스터를 형성하며, 특정 기능 단위로 작동하는 경향을 보임  
- **중심 유전자 식별**: PageRank, Degree 중심성 지표를 통해 네트워크 내에서 핵심적인 역할을 하는 유전자들을 확인함  


