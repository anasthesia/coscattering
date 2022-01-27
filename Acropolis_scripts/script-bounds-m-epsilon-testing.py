from acropolis.moreModels import ourDecayModel
import numpy as np

def main():
    deltaM=0.1
    logTanT=-4
    res=np.empty((0,29))
    for logEps in np.arange(-4.3,-4.2,0.2):
        for MD1 in np.logspace(1.7,2.5,num=7,base=10):
            print("Starting caclulation")
            #ourDecayModel(MD1,MD2,MZp,tan_theta,epsilon,alpha_X)
            model = ourDecayModel(MD1,MD1*(1+deltaM),3*MD1,10.**logTanT,10.**logEps,1/(4*np.pi))
            if model._tau() < 1e3:
                continue
            Yf = model.run_disintegration()
            res = np.append(res,np.array([np.concatenate(([MD1, logEps], Yf.flatten()))]),axis=0)
            print("Ending calculation")
    f=open("deltaM0.1_logTanT-4-testing-3.dat","w")
    np.savetxt(f,res)
    f.close()


if __name__ == "__main__":
    main()
