PriorityPruner
==============

Usage
-----
java -jar PriorityPruner.jar [OPTIONS]

Documentation
-------------
 * Documentation is available at [prioritypruner.sourceforge.net/documentation.html][1].


Building from source code
-------------------------

1) Install Maven v3 or higher
2) Install JDK v8 or higher
3) Extract source code
4) Install jar dependencies from source code folder

mvn install:install-file -Dfile=src/main/resources/lib/commons-cli-1.2.jar -DgroupId=commons-cli -DartifactId=commons-cli -Dversion=1.2 -Dpackaging=jar
mvn install:install-file -Dfile=src/main/resources/lib/log4j-1.2.17.jar -DartifactId=log4j -Dversion=1.2.17 -Dpackaging=jar

5) To run tests...

6) To build jar...


Developed By
============
 * Christopher K. Edlund
 * Malin Anker
 * Fredrick R. Schumacher
 * W. James Gauderman
 * David V. Conti
 

[1]: prioritypruner.sourceforge.net/documentation.html