import os
os.environ["PATH"] += os.pathsep + 'C:/Program Files/Graphviz/bin/'
from diagrams import Diagram
from diagrams import Cluster
from diagrams.aws.compute import *
from diagrams.aws.database import *
from diagrams.aws.network import ELB, Route53
from diagrams.aws.analytics import *
from diagrams.aws.business import *

with Diagram("Particle Collision", show=False) as web: # Keep
    collision = ComputeOptimizer("collision") # Keep
    LundString = Compute("Lund String") # Keep
    vphoton = Workmail("virtual photon") # Keep
    with Cluster("Initial Particles"): # Keep
        proton = A4B("proton")  # Keep
        electron = Chime("electron") # Keep
    with Cluster("Proton"): # Keep
    collision = ComputeOptimizer("collision")
        quark = AMI("Up")

        diquark = Compute("diquark 2212")

    Lund = Compute("Lund")
    with Cluster("Final State"):
        hadron0 = Compute("Pi+")

        hadron1 = Compute("Pi-")

        hadron2 = Compute("Proton")

    with Cluster("Not Measured"):
    proton >> diquark >> Lund
    proton >> quark >> collision >> Lund
web