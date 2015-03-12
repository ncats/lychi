LyChI (Layered Chemical Identifier)
===================================

This directory contains code for BARD structure standardizer. To build
simply type

```
mvn clean package
```

This produces two jar files: one with dependencies included, and without. The 
self-contained jar file can be invoked directly;
e.g., 

```
java -jar target/lychi-0.1-jar-with-dependencies.jar tests/standardizer_case1.smi
```

