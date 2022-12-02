import os
os.environ["PATH"] += os.pathsep + 'C:/Program Files/Graphviz/bin/'
from diagrams import Diagram
from diagrams import Cluster
from diagrams.aws.compute import *
from diagrams.aws.database import *
from diagrams.aws.network import ELB, Route53
from diagrams.aws.analytics import *
from diagrams.aws.business import *

with Diagram("Particle Collision Test", show=False): # Keep
    collision = ComputeOptimizer("collision") # Keep
    vphoton = Workmail("virtual photon") # Keep
    with Cluster("Initial Particles"): # Keep
        proton = A4B("proton")  # Keep
        electron = Chime("electron") # Keep
    with Cluster("Proton"): # Keep
        diquark = Compute("diquark 2101")

        quark = AMI("Up")

    Lund = Compute("Lund")
    struck_quark = AMI("Struck Up")

    with Cluster("Final State"):
        with Cluster("Hadrons"):
            hadron0 = Compute("Proton")

            hadron1 = Compute("Anti-Lambda")

            hadron2 = Compute("Sigma+")

            hadron3 = Compute("Pi-")

            hadron4 = Compute("Pi0")

            hadron5 = Compute("Pi+")

        with Cluster("Decay Photon Pairs"):
            fphoton0 = Neptune("Photon")

            fphoton1 = Neptune("Photon")

        scattered_electron = Chime("Scattered Electron")
    hadron4 >> fphoton1
    fphoton1 - fphoton0
    hadron1 >> hadron4
    hadron2 >> hadron5
    proton >> diquark >> Lund
    proton >> quark >> collision >> struck_quark >> Lund 
    Lund >> hadron0
    Lund >> hadron1
    Lund >> hadron2
    Lund >> hadron3
    electron >> vphoton >> collision
    electron >> scattered_electron
