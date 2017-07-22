#!/bin/bash
git pull origin master
git checkout a7784c2236544c55fdf4412f7b0f0719381215c6
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

cloneAndCheckout FangOost 09419c1388b70ecf69753f3e37f6276c3653e56f

cloneAndCheckout RungeKutta 6326974b245199852a1f00fccb4b677180ffc6d4

