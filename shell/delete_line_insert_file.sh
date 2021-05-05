for folder in disp-*; do cd ./$folder/; sed -i  '9,120d' 'POSCAR'; cd - ; done
for folder in ./disp-*; do cd ./$folder/; sed -i -e '8r ../Negorganic.vasp' 'POSCAR'; cd - ; done
