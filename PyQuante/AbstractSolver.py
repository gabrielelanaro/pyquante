from PyQuante import logging

class AbstractSolver:
    def __init__(self,molecule,**opts):
        logging.basicConfig(level=logging.INFO,format="%(message)s")
        self.molecule = molecule
        self.setup(**opts)
        self.print_setup_info()
        return

    def iterate(self,**opts):
        self.setup_iterations(**opts)
        self.print_pre_iteration_info()
        for self.iter in range(1,self.max_iter+1):
            self.iteration()
            if self.is_converged(): break
        self.print_post_iteration_info()
        return

    def setup_iterations(self,**opts):
        self.etol = opts.get('etol',1e-4)
        self.etemp = opts.get('etemp',None)
        self.max_iter = opts.get('max_iter',20)
        self.energy_history = []
        return

    def iteration(self):
        self.update_fock()
        self.solve_fock()
        self.calculate_energy()
        self.print_iteration_info()
        return

    def setup(self,**opts):
        self.method = "Abstract"
        self.nel = self.molecule.get_nel()
        self.nclosed,self.nopen = self.molecule.get_closedopen()
        self.Enuke = self.molecule.get_enuke()
        self.energy = 0
        return

    def print_setup_info(self):
        logging.info("\n%s calculation on %s" % (self.method,self.molecule.name))
        logging.debug("Nclosed = %d" % self.nclosed)
        logging.debug("Nopen = %d" % self.nopen)
        return

    def print_pre_iteration_info(self):
        logging.debug("Beginning SCF Optimization")
        logging.debug("It   Etot       Eone        Ej        Exc         Enuke")
        logging.debug("--   -------    -------     ------    -------     ------")
        return
    
    def print_post_iteration_info(self):
        if self.iter == self.max_iter-1:
            logging.warning("Warning! Maximum iterations (%d) reached" % self.max_iter)
        logging.info("%s final energy is %.6f reached in %d iterations"
                     % (self.molecule.name,self.energy,self.iter))
        return

    def met_target(self,etarget,targettol=1e-4):
        worked = abs(self.energy-etarget) < targettol
        if worked:
            logging.info("Target energy met")
        else:
            logging.warning("Warning! energy should have been %f" % etarget)
        return worked

    # These are effectively warnings to print out when something isn't
    # overloaded that should be.  Standard practice suggests that one
    # should actually raise exceptions here. I prefer this.
    def update_fock(self): print "Abstract: update fock"
    def solve_fock(self): print "Abstract: solve fock"
    def calculate_energy(self): print "Abstract: calculate energy"
    def is_converged(self):
        print "Abstract: is_converged"
        return True
    def print_iteration_info(self): print "Abstract: print iteration info"

