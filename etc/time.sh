cat config.py | grep multi_process > test.txt;
git show -q >> test.txt;
{ time python3 main.py >> test.txt & } 2>> test.txt; 
