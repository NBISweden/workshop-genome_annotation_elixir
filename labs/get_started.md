---
layout: default-overview
title: Get started
exercises: 10 min
questions:
  - How do I connect to server/VM?
  - How do I get the data and tools necessary for the course?
objectives:
  - Connect to server/VM
  - Get all the data necessary and the tools running for the course
---

## Foreword:

We will for all exercises use data for the fruit fly, *Drosophila melanogaster*, as that is one of the currently best annotated organisms and there is plenty of high quality data available. However, working on eukaryotes can be time consuming. Even a small genome like Drosophila would take too long to run within the time we have for this course. Thus to be sure to perform the practicals in good conditions, we will use the smallest chromosome of the drosophila (chromosome 4) like it was a whole genome.

## Prerequisites

### Connection to the server/VM 
Please connect yourself to the server/VM.

### Loading necessary tools  
All the tools are installed and you can load them by activating conda environment.

  ```bash
  conda activate gaas #or agat,....
  ```

   Now the GAAS environment is displayed at the beginnig of each prompt line: `(gaas)`
   To get out of the environment and restore your previous environment type:

  ```bash
  deactivate
  ```

### Setup your general working place    
Then you will prepare your general working place.  

   * Create your private place where all the magic will happen.  

   ```bash
   mkdir $USER
   cd $USER
   ```

   * get the data  
   Data needed for the exercices have to be copied locally.  

   ```bash
   cp -r /home/data/data_annotation/ .
   ```
