cat config.py | grep multi_process > main_output.txt;
git show -q >> main_output.txt;
{ time python3 main.py >> main_output.txt & } 2>> main_output.txt; 
