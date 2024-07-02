#!/bin/bash

# Set the classpath to include the current directory and the lib directory
CLASSPATH=.:./lib/*

# Compile the Java source file
javac -cp $CLASSPATH src/toroidalDiffusion/PDBPhiPsiExtractor.java -d .

# Run the compiled Java program
java -cp $CLASSPATH toroidalDiffusion.PDBPhiPsiExtractor $1 $2
