git pull origin master
git checkout b5612a768bd391973f9b60e29946775e54892da6



CALL :cloneAndCheckout "FunctionalUtilities" "cb634c343585267d51409a4af1871dde613f8b9c"
CALL :cloneAndCheckout "CharacteristicFunctions" "59d5f2e1789bdcfc8d861645ef3d6392de79832f"
CALL :cloneAndCheckout "FangOost" "d74b4e541aa5569ccd99192fe3d256e7d23b3883"
CALL :cloneAndCheckout "RungeKutta" "b286ce51f5e6957a59a5a49a3d1c60abdd765af4"
EXIT /B %ERRORLEVEL%

:cloneAndCheckout 
cd ..
IF EXIST %1  (
    cd %1
    git pull origin master
) ELSE (
    git clone https://github.com/phillyfan1138/%1
	cd %1
)
git checkout %2
cd ..
cd FFTOptionPricing
EXIT /B 0
