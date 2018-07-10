cat config.py | grep multi_process > main_output.txt;
git show -q >> main_output.txt;
{ time python3 main.py >> test.txt & } 2>> main_output.txt; 
