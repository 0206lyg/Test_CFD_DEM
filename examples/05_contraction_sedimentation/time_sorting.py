import os

# 파일 이름 설정
input_filename = "log.pimpleLYJHFDIBFoam"
output_filename = "extracted_time_data.txt"

# 현재 스크립트가 위치한 디렉토리 기준으로 경로 설정
current_dir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in locals() else os.getcwd()
input_path = os.path.join(current_dir, input_filename)
output_path = os.path.join(current_dir, output_filename)

try:
    with open(input_path, 'r', encoding='utf-8', errors='ignore') as infile, \
         open(output_path, 'w', encoding='utf-8') as outfile:
        
        for line in infile:
            stripped_line = line.strip()
            
            # 1. deltaT 라인인 경우 그대로 기록
            if stripped_line.startswith("deltaT ="):
                outfile.write(line)
                
            # 2. Time 라인인 경우 기록한 뒤, 한 블록이 끝났으므로 빈 줄(\n) 추가
            elif stripped_line.startswith("Time ="):
                outfile.write(line)
                outfile.write("\n")  # 다른 블록과의 구분을 위한 공백 라인 삽입
                
    print(f"추출 완료! 블록별로 구분되어 '{output_filename}' 파일에 저장되었습니다.")

except FileNotFoundError:
    print(f"오류: 같은 디렉토리 내에서 '{input_filename}' 파일을 찾을 수 없습니다.")