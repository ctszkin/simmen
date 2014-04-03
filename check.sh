Rscript roxygen.r
sleep 1
R CMD check simmen --as-cran
sleep 1
rm -r simmen/src-i386
sleep 1
rm -r simmen/src-x64
sleep 1
