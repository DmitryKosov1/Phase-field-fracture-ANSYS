# Phase-field-fracture-ANSYS
### The phase field model for fracture builds upon the pioneering thermo-dynamic framework established by Griffith, where crack growth will take place if a critical energy release rate is attained. We provide an efficient and robust implementation of the phase field method in the commercial finite element package ANSYS, enabling to model interactions and branching of cracks of arbitrary topological complexity.

### in the source files there is a library "USERELEMLIB.dll", to which you can connect by adding the ANS_USER_PATH variable to the environment variables, where the value is the path to the "USERELEMLIB.dll" file. The second file "APDL- 2D crack example.txt" is an APDL script describing the model parameters. Using this script you can repeat the task presented in the video below. To obtain the phase field, you must enter the command in ANSYS: "plnsol,curr". The third file "USERELEMLIB.f" is executable program code

https://github.com/DmitryKosov1/Phase-field-fracture-ANSYS/assets/131136458/fb983067-12bb-4860-b1f3-759e900211cb



