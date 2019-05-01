LyChI (Layered Chemical Identifier)
===================================

This directory contains code for BARD/GSRS structure standardizer and hash-code generator. 

To build simply type:

```
bash make.sh
```

You will need to have maven installed. Inside the `make.sh` script it simply adds the dependencies and calls `mvn clean package`.


This will build 2 jar files (one with dependencies and one without). They can be located in the `target` directory, and will look like:

```
target/lychi-0.5.1.jar  
target/lychi-0.5.1-jar-with-dependencies.jar
```

The self-contained jar file can be invoked directly. For example:

```
java -jar target/lychi-0.5.1-jar-with-dependencies.jar tests/standardizer_case1.smi
```


