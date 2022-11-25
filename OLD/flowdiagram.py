import os
os.environ["PATH"] += os.pathsep + 'C:/Program Files/Graphviz/bin/'
from diagrams import Diagram
from diagrams import Cluster
from diagrams.aws.compute import *
from diagrams.aws.database import *
from diagrams.aws.network import ELB, Route53
from diagrams.aws.analytics import *
from diagrams.aws.business import *

with Diagram("Particle Collision Test", show=False) as web: # Keep
    collision = ComputeOptimizer("collision") # Keep
    vphoton = Workmail("virtual photon") # Keep
    with Cluster("Initial Particles"): # Keep
        proton = A4B("proton")  # Keep
        electron = Chime("electron") # Keep
    with Cluster("Proton"): # Keep
        quark = AMI("Anti-Up")

        pairquark = AMI("Up")

        baryon = Compute("Proton")

    struck_quark = AMI("Struck Anti-Up")

    Lund = Compute("Lund")
    with Cluster("Final State"):
        with Cluster("Hadrons"):
            hadron0 = Compute("Rho0")

            hadron1 = Compute("Rho+")

            hadron2 = Compute("Pi-")

            hadron3 = Compute("Pi+")

            hadron4 = Compute("Pi-")

            hadron5 = Compute("Pi+")

            hadron6 = Compute("Pi0")

        with Cluster("Decay Photon Pairs"):
            fphoton0 = Neptune("Photon")

            fphoton1 = Neptune("Photon")

        scattered_electron = Chime("Scattered Electron")
    hadron6 >> fphoton1
    fphoton1 - fphoton0
    hadron0 >> hadron3
    hadron0 >> hadron4
    hadron1 >> hadron5
    hadron1 >> hadron6
    proton >> quark >> collision >> struck_quark >> Lund 
    proton >> pairquark >> Lund
    proton >> baryon >> Lund
    Lund >> hadron0
    Lund >> hadron1
    electron >> vphoton >> collision
    electron >> scattered_electron
web