This github repository is made by Emma Kok. 
The files are used for the bachelor thesis project.

Calculating_velocity.ipynb and Distance_travelled.py are made by Emma Kok.
The first file contains the code for calculating the velocity of the different type of fibers.
The second file contains the code for calculating the travelled distance of the different types of fibers.

The rest of the files are originally made by Proff. Matthew MacLeod, Dr. Antonia Praetorius and Dr. Maria del Prado Domercq.
The .txt files are the input files and are adjusted by Emma Kok. The rest of the files contains information and are also partially adjusted by Emma Kok.

More information about the model can be found underneath.

# Full-Multi-Microplastics-River-Model
=======================================

Microplastics River Model from Nano2Plast project

This repository contains one of the models developed in the ECO48-Nano2plast Project: EXTENDING NANOPARTICLE MODELS TO OPEN SOURCE MODELS OF THE FATE AND TRANSPORT OF MICROPLASTIC IN AQUATIC SYSTEMS.

The Full Multi Team: Proff. Matthew MacLeod (@MacLeodMatt), Dr. Antonia Praetorius (@apraetorius) and Dr. Maria del Prado Domercq (@PradoDomercq)

The code contained in this repository is intended as a framework to study microplastics fate and transport along a generic river system.

The model can be parameterised for different microplastic types (composition, shape, size, etc.) and river characteristics (flow velocity, suspended particulates, dimensions, etc.). Furthermore, this framework is intended as a flexible tool where new or updated fate processes can be included or reparameterised easily. 

This model follows a modular multimedia mass-balance modelling approach in which a system of coupled mass balance equations is built describing the transformation and transport processes of the plastic particles within the modelled system with specific rate constants. All processes are described using first-order kinetics based on our own review of the literature.

The modelled system consists of a generic one directional river structure where the river is subdivided in a set of river sections (RS) connected horizontally. Each river section is, at the same time, subdivided into four compartments of different types and dimensions representing a surface water layer (w0), a flowing water body (w1), a stagnant water body compartment (w2) and the top sediment layer of the river bed (sed). The system of differential equations is parameterized according to each compartment properties and dimensions and solved in a dynamic mode so that it yields the particles’s concentration (in number or mass per volume) in each river section and compartment as a function of time.

### Author of the model
===========
Maria del Prado Domercq
prado.domercq@gmail.com

### Adjuster of the model
Emma Kok
emma.kok@student.uva.nl
