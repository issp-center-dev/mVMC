x=1
end=10
while [ $x -lt $end ]
do
    cp ./StdFace_org.def ./StdFace.def
    echo "RndSeed = $x" >> StdFace.def
    ./vmc.out -s ./StdFace.def
    x=`expr $x + 12345`
    mv ./output ./output$x
done

cp ./StdFace_org.def ./StdFace_org.def.bak
rm ./*.def
rm ./*.dat
rm ./*.gp
mv ./StdFace_org.def.bak ./StdFace_org.def

