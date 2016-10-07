# -*- coding: utf-8 -*-


class RelaxationTensor:


    def secularize(self):
        """Secularizes the relaxation tensor


        """
        if self.as_operators:
            raise Exception("Cannot be secularized in an opeator form")
            
        else:
            N = self.data.shape[0]
            for ii in range(N):
                for jj in range(N):
                    for kk in range(N):
                        for ll in range(N):
                            if not (((ii == jj) and (kk == ll)) 
                                or ((ii == kk) and (jj == ll))) :
                                    self.data[ii,jj,kk,ll] = 0
                                        