PriorityPruner
==============

Usage
-----
`java -jar PriorityPruner.jar [OPTIONS]`


Documentation
-------------
Documentation is available at [chrisedlund.github.io/prioritypruner/documentation.html](chrisedlund.github.io/prioritypruner/documentation.html)


Building from source code
-------------------------

1) Install [Git](https://git-scm.com/)

2) Clone PriorityPruner code repository:

`git clone https://github.com/chrisedlund/prioritypruner.git prioritypruner`

3) Install [JDK](http://www.oracle.com/technetwork/java/javase/downloads/index.html) version 8 or higher, and make sure `$JAVA_HOME` points to this installation

4) Install [Maven](https://maven.apache.org/) version 3 or higher, and add the Maven bin folder to `$PATH`

5) Install jar dependencies from source code folder into your local Maven repository:

```
cd prioritypruner-code

mvn install:install-file -Dfile=src/main/resources/lib/commons-cli-1.2.jar \
   -DgroupId=commons-cli -DartifactId=commons-cli -Dversion=1.2 -Dpackaging=jar
   
mvn install:install-file -Dfile=src/main/resources/lib/log4j-1.2.17.jar \
   -DgroupId=log4j -DartifactId=log4j -Dversion=1.2.17 -Dpackaging=jar
```

6) You should now be ready to run various phases of the Maven lifecycle:

- To validate the project is correct and all necessary information is available, run:

`mvn validate`

- To compile the source code, run:

`mvn compile`

- To clean the build (delete everything in the target folder), run:

`mvn clean`

- To run unit and regressions tests, run:

`mvn test`

- To package as a jar (in the "target" folder), run:

`mvn package`


- To package as a jar (in the "target" folder) without running tests, run:

`mvn package -Dmaven.test.skip=true`

Developed By
============
 * Christopher K. Edlund [cedlund@usc.edu](mailto:cedlund@usc.edu)
 * Malin Anker
 * Fredrick R. Schumacher
 * W. James Gauderman
 * David V. Conti
 

