cat config.py | grep multi_process > out_main.txt;
git show -q >> out_main.txt;
{ time python3 main.py >> out_main.txt & } 2>> out_main.txt; 
