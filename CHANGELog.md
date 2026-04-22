이 문서는 각 한 주 단위로 개발 내용을 정리하기 위한 것이다.

# 2022년 10월 16일
issue #6
jquery 의 autocomplete 를 사용하게 테스트 코드를 작성하였다.
issue #16 도 마무리 하였음

# 2022년 10월 8일
issue #16
10월 4일 미팅때 논의 하였던 blast 에서 찾은 HSP 들에 해당하는 gene sequence 를 내려받기 위해서 Entrez 에서 gb 파일을 내려받고 파싱하는 부분을 구현.
결과물은 seqeuence 전체와 각 gene 에 대한 start, end, strand 를 파싱하는 것까지 구현하였다.
관련 코드는 `download_genebank.py` 에 있으며 테스트코드는 `tests/test_download_genebank.py` 에 있다.
테스트용 gb 파일은 `tests/blast_example/blast_example-CP054795.1.gb`에 저장하였다.
이 파일의 경로는 config.yaml 의 TEMPDIR 과 Common.get_gb_fn 에서 관리한다.

# 2022년 10월 3일 ~ 10월 9일
[마일스톤](https://github.com/io2oi/ctlabs-primerwrapper/milestone/1)
issue #6 작업 : 우선은 issue 에 적어 놓은데로 fullleneage 파일을 js 파일로 파싱해 놓아야 함.

# 2022년 9월 26일 ~ 2022년 10월 3일
issue #11, issue #12
git branch: ft-12-msa
multiple sequence alignment 에 대해서 mockup 을 만들었음.
로컬 컴퓨터에서 jalview 를 사용한다는 가정하에 만듦
