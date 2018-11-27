################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../CodonTable.cpp \
../CovarianceMatrix.cpp \
../FONSEModel.cpp \
../FONSEParameter.cpp \
../Gene.cpp \
../Genome.cpp \
../MCMCAlgorithm.cpp \
../Model.cpp \
../PAModel.cpp \
../PANSEModel.cpp \
../PANSEParameter.cpp \
../PAParameter.cpp \
../Parameter.cpp \
../RCPP_Model.cpp \
../RCPP_Parameter.cpp \
../RCPP_Trace.cpp \
../ROCModel.cpp \
../ROCParameter.cpp \
../SequenceSummary.cpp \
../Testing.cpp \
../Trace.cpp \
../main.cpp 

C_SRCS += \
../registerDynamicSymbol.c 

OBJS += \
./CodonTable.o \
./CovarianceMatrix.o \
./FONSEModel.o \
./FONSEParameter.o \
./Gene.o \
./Genome.o \
./MCMCAlgorithm.o \
./Model.o \
./PAModel.o \
./PANSEModel.o \
./PANSEParameter.o \
./PAParameter.o \
./Parameter.o \
./RCPP_Model.o \
./RCPP_Parameter.o \
./RCPP_Trace.o \
./ROCModel.o \
./ROCParameter.o \
./SequenceSummary.o \
./Testing.o \
./Trace.o \
./main.o \
./registerDynamicSymbol.o 

CPP_DEPS += \
./CodonTable.d \
./CovarianceMatrix.d \
./FONSEModel.d \
./FONSEParameter.d \
./Gene.d \
./Genome.d \
./MCMCAlgorithm.d \
./Model.d \
./PAModel.d \
./PANSEModel.d \
./PANSEParameter.d \
./PAParameter.d \
./Parameter.d \
./RCPP_Model.d \
./RCPP_Parameter.d \
./RCPP_Trace.d \
./ROCModel.d \
./ROCParameter.d \
./SequenceSummary.d \
./Testing.d \
./Trace.d \
./main.d 

C_DEPS += \
./registerDynamicSymbol.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


