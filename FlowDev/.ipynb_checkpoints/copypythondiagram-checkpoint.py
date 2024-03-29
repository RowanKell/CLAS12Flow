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

            hadron1 = Compute("Pi+")

            hadron2 = Compute("Pi-")

            hadron3 = Compute("test")

        scattered_electron = Chime("Scattered Electron")
    proton >> diquark >> Lund
    proton >> quark >> collision >> struck_quark >> Lund 
    electron >> vphoton >> collision
    electron >> scattered_electron
