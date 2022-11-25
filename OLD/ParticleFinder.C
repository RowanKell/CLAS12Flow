int ParticleFinder() {
    gROOT->ProcessLine("#include <vector>");
    
    auto hipofile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3301_3.hipo";
    HipoChain chain;
    
    chain.Add(hipofile);
    auto config_c12 = chain.GetC12Reader();
    
    config_c12->addExactPid(11,1);    //exactly 1 electron
    config_c12->addExactPid(111,1);    //exactly 1 pi0
    config_c12->addExactPid(2212,1);    //exactly 1 proton
    auto& c12 = chain.C12ref();
    
    int event_count = 0;
    // Loop over all events in Hipo file
    while (chain.Next()==true) {
        if (c12->getDetParticles().empty())
            continue;
        event_count += 1;
        if(event_count >= 30) {
            break;
        }
        event_count += 1;
        std::cout << "event #" << event_count << endl;
    }
    return 0;
}