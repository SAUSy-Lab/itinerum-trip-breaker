cat config.py | grep multi_process > ./etc/out_main.txt;
git show -q >> ./etc/out_main.txt;
{ time python3 main.py >> ./etc/out_main.txt & } 2>> ./etc/out_main.txt; 
