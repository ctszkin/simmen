R CMD build simmen
R CMD INSTALL simmen --build
sleep 1
rm -r simmen/src-i386
sleep 1
rm -r simmen/src-x64
sleep 1
