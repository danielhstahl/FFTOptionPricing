#!/bin/bash
git pull origin master
git checkout 2531e1bbbcf30889aacb82377868b02520f92da9
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

cloneAndCheckout FunctionalUtilities 33da1e6317022d0ccad5bb116bab6e1d4784439e

cloneAndCheckout CharacteristicFunctions 0dfefb3dc845d40ae94b70edb23cf6301bd6df18

cloneAndCheckout FangOost 3614c1e8a6f0608ba280ba25702cb6fc498e0428

cloneAndCheckout RungeKutta 6326974b245199852a1f00fccb4b677180ffc6d4

