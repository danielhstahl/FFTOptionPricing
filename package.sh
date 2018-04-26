#!/bin/bash
git pull origin master
git checkout b346f6968b36d0837e9b3bb1740b5d7a1d4854b3
function cloneAndCheckout {
	cd ..
	if [ -d "$1" ]; then
		cd $1
		git pull origin master
	else
		git clone https://github.com/phillyfan1138/$1
		cd $1
	fi

	git checkout $2
	cd ..
	cd FFTOptionPricing
}

cloneAndCheckout FunctionalUtilities 20bb2555ccb00539496860bf598eb466d348c043

cloneAndCheckout CharacteristicFunctions 0dfefb3dc845d40ae94b70edb23cf6301bd6df18

cloneAndCheckout FangOost 5919d52d5a59b317878b61115acd6cda3a91b97f

cloneAndCheckout RungeKutta 6326974b245199852a1f00fccb4b677180ffc6d4
cloneAndCheckout GaussNewton
cloneAndCheckout TupleUtilities
cloneAndCheckout AutoDiff

cd ..
git clone https://github.com/kthohr/optim ./optim
# build and install
cd ./optim
./configure -p -d
make
cd ..
cd FFTOptionPricing
