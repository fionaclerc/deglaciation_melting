#!/bin/sh

echo "Enter directory name:"
read dir

cd $dir

cp ../write*py .

pvpython write_pvd2txt_comp.py
pvpython write_pvd2txt_temp.py
pvpython write_pvd2txt_vel.py

mkdir out_files
cp out*txt out_files

