v=3.3.9

if [ ! -f eigen-$v.tar.bz2 ] ; then
  wget https://gitlab.com/libeigen/eigen/-/archive/$v/eigen-$v.tar.bz2
fi
rm -f Eigen
tar -xvjf eigen-$v.tar.bz2
ln -s eigen-$v Eigen
