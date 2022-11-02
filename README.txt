Pasos a seguir para correr el trabajo:

1. Escribir en la consola: make

2. Escribir en la consola ./tp1

3. Cuando el código ya se esté ejecutando deberá ingresar el nombre del archivo de texto del test que desee correr y luego el valor de p para dicho test.


En el main.c además se provee una función para correr todos los test de la catedra con la implementación de matrices ralas con CSR. Se debe llamar en el main a correr_test_catedra_CSR(). También hay una función para correr los test de la catedra y los test nuestros con ambas implementaciones de matrices ralas. Se debe llamar en el main a correr_tests_catedra_DOKS_CSR o correr_tests_nuestros_DOKS_CSR. 


Para correr los test nuestros se debe entrar al archivo matriz_rala_CSR.cpp y comentar la linea 8 y descomentar la 9. Esto se debe hacer tanto para matriz_rala_CSR.cpp como matriz_rala.cpp en la cual se tendrá que comentar la 6 y descomentar la 7. Esto es porque tuvimos problemas pasando el path completo por parametro. 

