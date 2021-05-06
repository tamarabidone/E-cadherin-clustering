package com.company;

import java.util.*;
import java.util.stream.Collectors;


/**
 * Here this module is to achieve the force interaction between these two layers of Cadherin.
 * @author yangchen
 */

public class Force_Interaction {
    // instance variables
    private ArrayList<Cadherin> cadList_1;
    private ArrayList<Cadherin> cadList_2;
    private ArrayList<Actin> actinList_1;
    private ArrayList<Actin> actinList_2;
    private double timestep;

    // constructor
    public Force_Interaction(ArrayList<Cadherin> cadList_1, ArrayList<Cadherin> cadList_2, ArrayList<Actin> actinList_1,
                             ArrayList<Actin> actinList_2, double timestep) {
        this.cadList_1 = cadList_1;
        this.cadList_2 = cadList_2;
        this.actinList_1 = actinList_1; // get the actinList
        this.actinList_2 = actinList_2;
        this.timestep = timestep;
    }

    private double distance(Cadherin cad1, Cadherin cad2) {
        return Math.sqrt(Math.pow(cad1.getBead_position()[0] - cad2.getBead_position()[0], 2) +
                Math.pow(cad1.getBead_position()[1] - cad2.getBead_position()[1], 2));
    }

    private double distance(Cadherin cad1, Actin actin1) {
        return Math.sqrt(Math.pow(cad1.getBead_position()[0] - actin1.getBead_position()[0], 2) +
                Math.pow(cad1.getBead_position()[1] - actin1.getBead_position()[1], 2));
    }

    private void boundForceUpdate(Cadherin cad1, Cadherin cad2) {

        double k_it = 2e3;
        double[] tmp = new double[2];
        double[] tmp2 = new double[2];
        for (int i = 0; i < 2; i++) {
            tmp[i] = cad1.getBead_position()[i] - k_it*(cad1.getBead_position()[i] - cad2.getBead_position()[i])*timestep/cad1.getFrictional();
            tmp2[i] = cad2.getBead_position()[i] - k_it*(cad2.getBead_position()[i] - cad1.getBead_position()[i])*timestep/cad2.getFrictional();
        }
        cad1.setBead_position_x(tmp[0]); cad1.setBead_position_y(tmp[1]);
        cad2.setBead_position_x(tmp2[0]); cad2.setBead_position_y(tmp2[1]);
    }

    private void boundForceUpdate(Cadherin cad, Actin actin) {

        double k_it = 2e3;
        double[] tmp = new double[2];
        double[] tmp2 = new double[2];

        for (int i = 0; i < 2; i++) {
            tmp[i] = cad.getBead_position()[i] - k_it*(cad.getBead_position()[i] - actin.getBead_position()[i])*timestep/cad.getFrictional();
            tmp2[i] = actin.getBead_position()[i] - k_it*(actin.getBead_position()[i] - cad.getBead_position()[i])*timestep/actin.getFrictional();
        }
        cad.setBead_position_x(tmp[0]); cad.setBead_position_y(tmp[1]);
        actin.setBead_position_x(tmp2[0]); actin.setBead_position_y(tmp2[1]);
    }


    private void backToOrigin(Cadherin cad1, Cadherin boundCad) {

        // back to original setting
        cad1.setState(0);
        boundCad.setState(0);
        cad1.setForce(0);
        boundCad.setForce(0);
        cad1.setCadBoundIndex(-1);
        boundCad.setCadBoundIndex(-1);
    }

    // the catch bound function for the X-dimer (unbound)
    private double xCatch (double force) {
        return 100*1.0*Math.exp(-force*0.34/4.114)+0.7*Math.exp(force*0.34/4.114);
    }


    private double sSlip(double force) {
        return 1.27e-4*Math.exp(0.52*force);
    }

    /**
     * The implementation of trans-interaction of cad
     */
    public void trans_interaction() {// later may be I need to change this to trans-interaction

        //todo: double-check if there is any problem of recording force in between!!

        double cutoff = 0.003; // unit: micrometers, use the 15 nm as the cutoff distance
        double cutoff_2 = 0.001; // unit: micrometers

        // all kinetic rate constant regarding the transition between different states
        double k_fx = 3.8e4; // m -> x unit: 1/s
        double k_bx = 1.84e3; // x -> m unit: 1/s
        double k_fms = 3.1e-1; // m -> s unit: 1/s
        double k_bms = 1.27e-4; // s -> m unit: 1/s
        double k_fxs = 86; // x -> s unit: 1/s
        double k_bxs = 0.8; //s -> x unit: 1/s

        // using k * dt to calculate the probability of transition between different states
        double p_fx = k_fx * timestep;
        double p_bx = k_bx * timestep; // ideal
        double p_fms = k_fms * timestep;
        double p_bms = k_bms * timestep;
        double p_fxs = k_fxs * timestep;
        double p_bxs = k_bxs * timestep;
        //double p_rotate = 0.0341; // the probability that the range cadherin is in the range of rotation binding (1-cos(degree), degree = 15)

        // the spring constant we use to calculate the force when binding state is 1 (X)
        double k_it = 2e3; // unit: 2pn/um (ref!!)

        // After the update of corresponding state thermal diffusion and domain-diffusion
        // Here, consider the kinetic interaction for different states
        for (Cadherin cad1 : cadList_1) {
            // case 1: if the current state of cad1 is 0. (1) can bind to X-dimer (2) S dimer (3) remain monomer
            if (cad1.getState() == 0) {
                // iterate through all unbound cad in the layer 2 check the certain criteria
                for (Cadherin cad2 : cadList_2) {
                    if (cad2.getState() == 0) {
                        double dist = distance(cad1, cad2); // calculate the distance

                        // bound with cadherin in the layer 2 to form X-dimer
                        // (a) dist less than cut-off criteria (b) probability satisfied
                        if ((dist <= cutoff) && (Math.random() <= p_fx)) {

                            cad1.setState(1);
                            cad1.setCadBoundIndex(cadList_2.indexOf(cad2)); // set the value of index of bounded cadherin
                            cad2.setState(1);
                            cad2.setCadBoundIndex(cadList_1.indexOf(cad1)); // actually there is no need for this line

                            // update the position of the cadherin
                            boundForceUpdate(cad1, cad2);
                            double ndist = distance(cad1, cad2);
                            cad1.setForce(k_it * ndist);
                            cad2.setForce(k_it * ndist);


                            // break the inner for loop (which indicate the form of binding between two layers)
                            break;
                        }

                        if ((dist <= cutoff_2) && (Math.random() <= p_fms)) {
                            // from the S-dimer directly
                            // at this stage, we set the two caderin at the same position (x,y) in two layers
                            double force = k_it * dist;

                            cad1.setState(2);
                            cad2.setState(2);
                            cad1.setForce(force);
                            cad2.setForce(force);
                            cad1.setCadBoundIndex(cadList_2.indexOf(cad2));
                            cad2.setCadBoundIndex(cadList_1.indexOf(cad1));

                            break;
                        }
                    }
                }
                continue; // make sure we don't continue the rest this makes sure each cadherin in a specific state would only be tested the state once
            }


            // case 2: if the current state of the cadherin is 1 then (1) change to S-dimer (2) back to M (3) remain X
            if (cad1.getState() == 1) {
                // indicate cad2 / trbound must be 1
                // recalculate the force between them
                Cadherin boundCad = cadList_2.get(cad1.getCadBoundIndex());
                double ndist = distance(cad1, boundCad);// get corresponding bound cad in layer 2
                double nforce = k_it * ndist;
                if (Math.random() > xCatch(nforce)*timestep) {
                    boundForceUpdate(cad1, boundCad); // update the position
                    double nndist = distance(cad1, boundCad);
                    double nnforce = k_it*nndist;
                    cad1.setForce(nnforce);
                    boundCad.setForce(nnforce);

                    if ((nnforce <= 2) && (Math.random() <= p_fxs)) { // from X -> S (maybe make the cutoff force bigger)
                        cad1.setState(2);
                        boundCad.setState(2);
                    }
                } else {
                    backToOrigin(cad1, boundCad);
                }

                continue;
            }

            // case 3 if the current state is 2 (a) move back to M (b) move back to X (3) remain
            if (cad1.getState() == 2) {
                // (a) move back to the M
                Cadherin boundCad = cadList_2.get(cad1.getCadBoundIndex());
                // Still we can test the cadboundIndex here
                if (Math.random() <= p_bms) {
                    backToOrigin(cad1, boundCad);
                }

                if (Math.random() <= p_bxs) { // move back to X
                    cad1.setState(1);
                    boundCad.setState(1);
                }
            }
        }
    }

    /**
     * This function is to test whether the monomer in the cis-cluster will form the X or S dimer
     * depending on the distance and the probability
     */
    public void trans_cis_interaction() {
        // Here we only need to consider the monomer in the cis cluster
        double cutoff = 0.003; // unit: micrometers, use the 15 nm as the cutoff distance
        double cutoff_2 = 0.001; // unit: micrometers
        double k_fx = 3.8e4;
        double k_fms = 3.1e-1; // m -> s unit: 1/s
        double p_fx = k_fx * timestep;
        double p_fms = k_fms * timestep;
        double k_it = 2e3;
        double actink_it = 2e3;

        // filter all monomers in the cadlist1 and cadlist2
        List<Cadherin> cisMList1 = cadList_1.stream().filter(cad -> (cad.getState() == 3 && cad.isCisM())).collect(Collectors.toList());
        List<Cadherin> cisMlist2 = cadList_2.stream().filter(cad -> (cad.getState() == 3 && cad.isCisM())).collect(Collectors.toList());
        for (Cadherin cad1 : cisMList1) {
            for (Cadherin cad2 : cisMlist2) {
                double dist = distance(cad1, cad2);
                // form the X-dimer in the cis cluster
                if (dist <= cutoff && Math.random() <= p_fx) {
                    // update the corresponding state
                    cad1.setCisM(false);
                    cad2.setCisM(false);
                    cad1.setCisX(true);
                    cad2.setCisX(true);
                    cad1.setCadBoundIndex(cadList_2.indexOf(cad2));
                    cad2.setCadBoundIndex(cadList_1.indexOf(cad1));
                    boundForceUpdate(cad1, cad2);
                    if (cad1.getActState() > 0) {
                        Actin actin1 = actinList_1.get(cad1.getActBoundIndex());
                        boundForceUpdate(cad1, actin1);
                        actin1.setForce(actink_it * distance(cad1, actin1));
                    }
                    if (cad2.getActState() > 0) {
                        Actin actin2 = actinList_2.get(cad2.getActBoundIndex());
                        boundForceUpdate(cad2, actin2);
                        actin2.setForce(actink_it * distance(cad2, actin2));
                    }
                    double ndist = distance(cad1, cad2);
                    cad1.setForce(k_it * ndist);
                    cad2.setForce(k_it * ndist);

                    break;
                }

                // form the S-dimer directly
                if (dist <= cutoff_2 && Math.random() <= p_fms) {
                    cad1.setCisM(false);
                    cad2.setCisM(false);
                    cad1.setCadBoundIndex(cadList_2.indexOf(cad2));
                    cad2.setCadBoundIndex(cadList_1.indexOf(cad1));
                    boundForceUpdate(cad1, cad2/*, dist*/);
                    if (cad1.getActState() > 0) {
                        Actin actin1 = actinList_1.get(cad1.getActBoundIndex());
                        boundForceUpdate(cad1, actin1/*, distance(cad1, actin1)*/);
                        actin1.setForce(actink_it * distance(cad1, actin1));
                    }
                    if (cad2.getActState() > 0) {
                        Actin actin2 = actinList_2.get(cad2.getActBoundIndex());
                        boundForceUpdate(cad2, actin2/*, distance(cad2, actin2)*/);
                        actin2.setForce(actink_it * distance(cad2, actin2));
                    }
                    double ndist = distance(cad1, cad2);
                    cad1.setForce(k_it * ndist);
                    cad2.setForce(k_it * ndist);

                    break;
                }
            }
        }
    }

    /**
     * Here, we implement the third version of cis_interaction to simplify the forward process cis
     */
    public void cis_interaction3() {

        double k_cis = 100; // from s-dimer to cis-interface (I use the largest one)
        double k_bcis = 0.1; // from cis-interface back to S-dimer

        // probability to achieve the transition
        double p_cis = k_cis * timestep;
        double p_bcis = k_bcis * timestep;
        double cutoff_cis = 0.008;  // unit: um from the paper to form the a cis-interface
        double probability = 0.001; // the probability to break the monomer in the cis interaction
        /**
         * Here we implement the cis backward
         * (1) decide to break the cis bond or not in the corresponding layer separately
         * (2) iterate through to check all 3.0 cases
         */
        // determine the first layer
        cisBackward(cadList_1, p_bcis, probability); // must have some (3,0) state cads
        // determine the seconde layer
        cisBackward(cadList_2, p_bcis, probability); // must have some (3,0) state cads

        // Step 2: we iterate through the cadList1 and 2 to update the 3.0 state cad
        // layer one as the ref layer
        for (Cadherin cad1 : cadList_1) {
            if (cad1.getState() == 3 && cad1.getCisState() == 0) {
                if (cad1.isCisM()) {
                    cad1.setState(0);
                    cad1.setCadBoundIndex(-1);
                    cad1.setForce(0);
                    cad1.setCisM(false);
                    cad1.setCisX(false);
                    if (cad1.getActState() > 0) {
                        actinList_1.set(cad1.getActBoundIndex(), null);
                        cad1.setActState(0);
                    }
                } else {
                    // it is a dimer
                    Cadherin trBound1 = cadList_2.get(cad1.getCadBoundIndex());
                    if (trBound1.getCisState() == 0) {
                        // actin will be removed
                        if (cad1.getActState() > 0) {
                            actinList_1.set(cad1.getActBoundIndex(), null);
                            cad1.setActState(0);
                        }

                        if (trBound1.getActState() > 0) {
                            actinList_2.set(trBound1.getActBoundIndex(), null);
                            trBound1.setActState(0);
                        }
                        if (cad1.isCisX()) {
                            // it is a X dimer
                            cad1.setState(1);
                            trBound1.setState(1);
                            cad1.setCisX(false);
                            trBound1.setCisX(false);
                        } else {
                            // it is a S dimer
                            cad1.setState(2);
                            trBound1.setState(2);
                        }
                    }
                }
            }
        }

        // step3: update the state for 3.0 Monomer in the layer two
        for (Cadherin cad2 : cadList_2) {
            if (cad2.getState() == 3 && cad2.getCisState() == 0 && cad2.isCisM()) {
                cad2.setState(0);
                cad2.setCadBoundIndex(-1);
                cad2.setForce(0);
                cad2.setCisM(false);
                cad2.setCisX(false);
                if (cad2.getActState() > 0) {
                    actinList_2.set(cad2.getActBoundIndex(), null);
                    cad2.setActState(0);
                }
            }
        }


        /**
         * Here, we implement the cisforward function to only consider the formation in same layer
         */
        // forward cis interaction in the layer one
        cisForward(cadList_1, p_cis, cutoff_cis);
        // forward cis interaction in the layer two
        cisForward(cadList_2, p_cis, cutoff_cis); // in the same layer




        // deal with the boundary case when only one side of 2.0 form the cis, we update the other side to 3.0
        List<Cadherin> cad2layer1 = cadList_1.stream().filter(cad -> cad.getState() == 2).collect(Collectors.toList());
        for (Cadherin cad : cad2layer1) {
            Cadherin trBound = cadList_2.get(cad.getCadBoundIndex());
            if (trBound.getState() == 3) {
                cad.setState(3);
            }
        }

        // same for layer2
        List<Cadherin> cad2layer2 = cadList_2.stream().filter(cad -> cad.getState() == 2).collect(Collectors.toList()); // here only two types of S (1) from trans_interaction (2) from cisS backward
        for (Cadherin cad : cad2layer2) {
            try {
                Cadherin trBound = cadList_1.get(cad.getCadBoundIndex());
                if (trBound.getState() == 3) {
                    cad.setState(3);
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    private void cisForward(ArrayList<Cadherin> cadList, double p_cis, double cutoff_cis) {

        // O(n^2)
        // here we implement the forward process of cis_interaction
        // step1: filter out the cadherin to form cis_interaction
        List<Cadherin> cisCads = cadList.stream().filter(cad -> (cad.getState() == 2 || (cad.getState() == 3 && cad.getCisState() < 2))).collect(Collectors.toList());
        int size = cisCads.size();
        for (int i = 0; i < size; i++) {
            Cadherin cad1 = cisCads.get(i); // note this i is not the cad position in the cadlist arraylist
            if (cad1.getCisState() < 2) {
                int index1 = cadList.indexOf(cad1);
                for (int j = i + 1; j < size; j++) {
                    Cadherin cad2 = cisCads.get(j);
                    if (cad2.getCisState() < 2) { // still can add new index
                        int index2 = cadList.indexOf(cad2);
                        if ((!cad1.getCadCisBoundIndex().contains(index2)) && (!cad2.getCadCisBoundIndex().contains(index1))) {
                            double distance = distance(cad1, cad2);
                            if ((distance <= cutoff_cis) && (Math.random() <= p_cis)) {
                                // then we say they can form this cis_interaction and update the state
                                //1. add the index
                                cad1.getCadCisBoundIndex().add(index2);
                                cad2.getCadCisBoundIndex().add(index1);
                                //2. change the cisstate
                                cad1.setCisState(cad1.getCadCisBoundIndex().size());
                                cad2.setCisState(cad2.getCadCisBoundIndex().size());
                                //3. update the state
                                cad1.setState(3);
                                cad2.setState(3);
                            }
                        }
                    }
                }
            }
        }
    }

    private void cisBackward(ArrayList<Cadherin> cadList, double p_bcis, double probability) {

        // in this method, we determine the break cis bond behaviour
        int size = cadList.size(); // this part write very complicated
        for (int i = 0; i < size; i++) {
            Cadherin cad = cadList.get(i);
            if (cad.getState() == 3 && cad.getCisState() > 0) { // with cis cad in this layer
                int index1 = cadList.indexOf(cad); // this is i actually
                int size2 = cad.getCadCisBoundIndex().size(); // how many cis cad bound
                for (int j = size2 - 1; j >= 0; j--) {
                    int index = cad.getCadCisBoundIndex().get(j);
                    if (!cad.isCisM()) { // cisX or cisS
                        if (index > index1) {
                            // the ciscad is behind the cad, therefore we haven't check this bond
                            Cadherin ciscad = cadList.get(index);
                            if (!ciscad.isCisM()) {
                                if (Math.random() <= p_bcis) {
                                    // remove the index
                                    ciscad.getCadCisBoundIndex().remove(Integer.valueOf(index1));
                                    cad.getCadCisBoundIndex().remove(Integer.valueOf(index));
                                    // update state
                                    ciscad.setCisState(ciscad.getCadCisBoundIndex().size());
                                }
                            } else {
                                // the ciscad is a monomer
                                if (Math.random() <= probability) {
                                    ciscad.getCadCisBoundIndex().remove(Integer.valueOf(index1));
                                    cad.getCadCisBoundIndex().remove(Integer.valueOf(index));
                                    // update the state of ciscad
                                    ciscad.setCisState(ciscad.getCadCisBoundIndex().size());
                                }
                            }
                        }
                    } else {
                        if (index > index1) {
                            Cadherin ciscad = cadList.get(index);
                            if (Math.random() <= probability) {
                                ciscad.getCadCisBoundIndex().remove(Integer.valueOf(index1));
                                cad.getCadCisBoundIndex().remove(Integer.valueOf(index));
                                // update the state of ciscad
                                ciscad.setCisState(ciscad.getCadCisBoundIndex().size());
                            }
                        }
                    }
                }
                // update the state of cad
                cad.setCisState(cad.getCadCisBoundIndex().size());
            }
        }
    }

    // a series of method of catch bond curves
    // from weak to strong
    private double f12(double force) {
        return 3*Math.exp(force*0.2 / 4.114);
    }

    // from strong to weak
    private double f21(double force) {
        return 20*Math.exp(-force*4 / 4.114);
    }

    // Strong catch bond curve
    private double fstrong(double force) {
        return Math.exp(-force*2/4.114) + 0.0003*Math.exp(force*2/4.114);
    }

    // Weak catch bond curve
    private double fweak(double force) {
        return 5*Math.exp(-force*2/4.114) + 0.003*Math.exp(force*2/4.114);
    }

    /**
     * The third implementation of actin interaction
     * This is a intermediate method between the first and the second method
     */
    public ArrayList<Map<Integer, Integer>> actin_interaction3() {

        // we wanna have the corresponding position
        double zAct_1 = -0.01;
        double zAct_2 = 0.03;
        int size = cadList_1.size();
        double k_it = 2e3;
        double actink_it = 2e3;
        double force = 50; // force to remove the corresponding cad in the same level
        double probability = 1e-3; // later we test different binding probability

        // the add process of actin will not be affected in the way we apply the actin
        // First, we consider how do we add Actin for the two corresponding cad layer
        // Because later these two layer will have different situation based on the actin movement
        // Then, we iterate through cad list independently!
        // This part is easy to modify, we just iteraction through and add the actin! doesn't care about the real situation
        for (int i = 0; i < size ; i++) { // we can still use the filter to do this step, but for now we use for loop
            Cadherin cad1 = cadList_1.get(i);
            Cadherin cad2 = cadList_2.get(i);

            // this step include the addtion of actin to each cad in the cis interaction each iteration according to the probability
            if (cad1.getState() == 3 && cad1.getActState() == 0 && (Math.random() <= probability)) {
                double[] tmp = cad1.getBead_position().clone();
                tmp[2] = zAct_1;
                Actin actin = new Actin(tmp);
                // record the corresponding index of cad to the actin
                actin.setBoundcadIndex(i);
                // notice that later when we delete the act just assign that position to null and no need to reset the actBoundindex = -1
                if (-1 == cad1.getActBoundIndex()) { // first-time add the actin
                    actinList_1.add(actin);
                    // record the corresponding index of actin to the cad
                    cad1.setActBoundIndex(actinList_1.indexOf(actin));
                } else { // later times replace the corresponding position in the actin Arraylist
                    actinList_1.set(cad1.getActBoundIndex(), actin);
                }
                cad1.setActState(1); // initially we think all bond are weakly bond
            }

            if (cad2.getState() == 3  && (Math.random() <= probability) && cad2.getActState() == 0) {
                double[] tmp = cad2.getBead_position().clone();
                tmp[2] = zAct_2;
                Actin actin = new Actin(tmp);
                // record the index of cad
                actin.setBoundcadIndex(i);
                if (-1 == cad2.getActBoundIndex()) {
                    actinList_2.add(actin);
                    // record the actin
                    cad2.setActBoundIndex(actinList_2.indexOf(actin));
                } else {
                    actinList_2.set(cad2.getActBoundIndex(), actin);
                }
                cad2.setActState(1);
            }
        }

        // Then we implement the new logic of finding the small clusters and move cads accordingly
        // move the first layer
        ArrayList<HashSet<Integer>> Clusters1 = moveCisCad2(cadList_1, actinList_1, force);
        // move the second layer
        ArrayList<HashSet<Integer>> Clusters2 = moveCisCad2(cadList_2, actinList_2, force);

        ArrayList<Map<Integer, Integer>> info = ClusterInfo(Clusters1, Clusters2);

        // To print out the corresponding two-layer cluster info, we don't do the boundary check
        Map<Integer, Integer> twolayer_info = updateBoundary(Clusters1, Clusters2); // we do the statistics even when dont need to filp it
        info.add(twolayer_info); // add the two layer info to the info arraylist

        // decide the fate
        // here we assume that if force = 0 or prop = 0 the actin does not affect the trans cadherin
        if (force != 0 && probability != 0) {
            // filter the X and S dimer in the cis
            List<Cadherin> cisCad = cadList_1.stream().filter(cad -> (cad.getState() == 3) && (!cad.isCisM())).collect(Collectors.toList());
            for (Cadherin cad : cisCad) {
                Cadherin trBound = cadList_2.get(cad.getCadBoundIndex());
                if (cad.getActBoundIndex() != -1 && trBound.getActBoundIndex() != -1) { // todo: need double check
                    deterFate(cad, trBound, actinList_1.get(cad.getActBoundIndex()), actinList_2.get(trBound.getActBoundIndex()),actink_it, k_it);
                }
            }
        } else {
            // deal with the case of no force and no binding, which means the actin does not affect the cis cluster
            // we simply update the relative position for cadherins in the cis cluster
            // then when we iterate through the first layer the second layer must exist
            for (Cadherin cad : cadList_1) {
                // iterate through the cadlist 1
                if (cad.getState() == 3) {
                    // update the position
                    Cadherin boundCad = cadList_2.get(cad.getCadBoundIndex());
                    boundForceUpdate(cad, boundCad);
                    cad.setForce(k_it*distance(cad, boundCad));
                    boundCad.setForce(k_it*distance(cad, boundCad));
                }
            }

        }

        return info;
    }

    private ArrayList<Map<Integer, Integer>> ClusterInfo(ArrayList<HashSet<Integer>> Clusters1, ArrayList<HashSet<Integer>> Clusters2) {

        // here we do some statistics to get the info of each Clusters
        ArrayList<Map<Integer, Integer>> info = new ArrayList<>(2);
        Map<Integer, Integer> info1 = new HashMap<>();
        Map<Integer, Integer> info2 = new HashMap<>();
        int size1 = (int) cadList_1.stream().filter(cad -> (cad.getState() == 3 && cad.getCisState() == 0)).count();
        int size2 = (int) cadList_2.stream().filter(cad -> (cad.getState() == 3 && cad.getCisState() == 0)).count();
        Count(Clusters1, info1, size1);
        Count(Clusters2, info2, size2);

        info.add(info1);
        info.add(info2);

        return info;
    }

    private void Count(ArrayList<HashSet<Integer>> Clusters, Map<Integer, Integer> info, int size) {

        // here we count count the info of cluster
        for (HashSet<Integer> cluster : Clusters) {
            int clusterSize = cluster.size();
            if (clusterSize != 0) {
                if (!info.containsKey(clusterSize)) {
                    info.put(clusterSize, 1);
                } else {
                    info.put(clusterSize, info.get(clusterSize) + 1);
                }
            }
        }
        info.put(1, size);
    }

    private Map<Integer, Integer> updateBoundary(ArrayList<HashSet<Integer>> clusters1, ArrayList<HashSet<Integer>> clusters2) {

        // in this method, we update the boundary of the cad in cis cluster if they cross the boundary as well as actin
        Map<Integer, Integer> twolayer_info = new HashMap<>();
        List<Cadherin> Singles1 = cadList_1.stream().filter(cad -> (cad.getState() == 3 && cad.getCisState() == 0)).collect(Collectors.toList());
        List<Cadherin> Singles2 = cadList_2.stream().filter(cad -> (cad.getState() == 3 && cad.getCisState() == 0)).collect(Collectors.toList());

        // we start from the layer1 singles
        findUpdate(Singles1, Singles2, clusters1, clusters2, 1, twolayer_info);
        // then we start from the layer 2 singles
        findUpdate(Singles2, Singles1, clusters1, clusters2, 2, twolayer_info);

        findUpdate2(clusters1, clusters2, twolayer_info);

        return twolayer_info;
    }

    private void findUpdate2(ArrayList<HashSet<Integer>> clusters1, ArrayList<HashSet<Integer>> clusters2, Map<Integer, Integer> info) {

        // this method we deal with all 2-layer cluster starting from the first layer
        double domainSize = 1;
        while (clusters1.size() > 1) {// to ensure we always get the last cad in the cluster
            // step 1: record the cluster + trans in two layers starting from 3.0 X, S in layer1
            ArrayList<Cadherin> crossCluster = new ArrayList<>(); // store the cadherin in this two-layer cluster
            ArrayList<Actin> actinCluster = new ArrayList<>(); // store the actins bind
            ArrayList<Integer> nextLevelCadindex = new ArrayList<>();
            ArrayList<Integer> curreentLevelCadindex = new ArrayList<>();
            boolean jump = true;
            boolean level = false;
            int crossSize = 0;
            int i = clusters1.size() - 1;
            for (int index : clusters1.get(i)) {
                Cadherin cad = cadList_1.get(index);
                crossCluster.add(cad);
                if (cad.getActState() > 0) actinCluster.add(actinList_1.get(cad.getActBoundIndex()));
                if (testBoundary(cad, domainSize)) crossSize++;

                if (!cad.isCisM()) {
                    int transindex = cad.getCadBoundIndex();
                    curreentLevelCadindex.add(transindex);
                }
            }
            clusters1.remove(clusters1.get(i)); // delete current cad from the singles

            // jump between the first and second layer
            while (jump) {
                nextLevelCadindex.clear();
                for (int index : curreentLevelCadindex) {
                    if (level) {
                        // here we search the first clusters
                        Cadherin trBound = cadList_1.get(index);
                        if (!crossCluster.contains(trBound)) {
                            crossCluster.add(trBound);
                            if (trBound.getActState() > 0) actinCluster.add(actinList_1.get(trBound.getActBoundIndex()));
                            if (testBoundary(trBound, domainSize)) crossSize++;
                            for (int j = clusters1.size() - 1; j > 0; j--) {
                                if (clusters1.get(j).contains(index)) {
                                    for (int cisindex : clusters1.get(j)) {
                                        if (cisindex != index) {
                                            Cadherin ciscad = cadList_1.get(cisindex);
                                            crossCluster.add(ciscad);
                                            if (ciscad.getActState() > 0) actinCluster.add(actinList_1.get(ciscad.getActBoundIndex()));
                                            if (testBoundary(ciscad, domainSize)) crossSize++;

                                            if (!ciscad.isCisM()) {
                                                int indexcad2 = ciscad.getCadBoundIndex();
                                                Cadherin cad2 = cadList_2.get(indexcad2);
                                                if (!crossCluster.contains(cad2)) {
                                                    nextLevelCadindex.add(indexcad2);
                                                }
                                            }
                                        }
                                    }
                                    // remove this cluster in the cluster2
                                    clusters1.remove(clusters1.get(j));
                                    break;
                                }
                            }
                        }
                    } else {
                        // search the second clusters
                        Cadherin trBound = cadList_2.get(index);
                        if (!crossCluster.contains(trBound)) {
                            // to avoid replication
                            crossCluster.add(trBound);
                            if (trBound.getActState() > 0) actinCluster.add(actinList_2.get(trBound.getActBoundIndex()));
                            if (testBoundary(trBound, domainSize)) crossSize++;
                            for (int j = clusters2.size() - 1; j > 0; j--) {
                                if (clusters2.get(j).contains(index)) {
                                    for (int cisindex : clusters2.get(j)) {
                                        if (cisindex != index) {
                                            Cadherin ciscad = cadList_2.get(cisindex);
                                            crossCluster.add(ciscad);
                                            if (ciscad.getActState() > 0) actinCluster.add(actinList_2.get(ciscad.getActBoundIndex()));
                                            if (testBoundary(ciscad, domainSize)) crossSize++;

                                            if (!ciscad.isCisM()) {
                                                int indexcad2 = ciscad.getCadBoundIndex();
                                                Cadherin cad2 = cadList_1.get(indexcad2);
                                                if (!crossCluster.contains(cad2)) {
                                                    // if the corresponding cad2 not contain in the crossCluster, indicating that it must be 3.1
                                                    nextLevelCadindex.add(indexcad2);
                                                }
                                            }
                                        }
                                    }
                                    // remove this cluster in the cluster2
                                    clusters2.remove(clusters2.get(j));
                                    break;
                                }
                            }
                        }
                    }
                }

                // here we iterate through the cluster contains the index in the currentlevel
                // we update several condtions
                if (nextLevelCadindex.size() == 0) {
                    jump = false;
                } else {
                    curreentLevelCadindex = (ArrayList<Integer>)(nextLevelCadindex.clone());
                    level = !level; // jump to the next level
                }
            }

            // count the two-layer cluster size and put into the map
            int twolayer_size = crossCluster.size();
            if (!info.containsKey(twolayer_size)) {
                info.put(twolayer_size, 1);
            } else {
                info.put(twolayer_size, info.get(twolayer_size) + 1);
            }

            // here we find the corresponding two level cluster
            // decide flip or not
            // we only move the cad and actin
            if (crossSize > 0) {
                if (crossSize == crossCluster.size()) {
                    // this indicates all cad in this cluster out of the boundary
                    // under this situation, we just do the normal boundary update
                    for (Cadherin cad : crossCluster) {
                        cad.boundary_force_update(domainSize);
                    }
                    for (Actin actin : actinCluster) {
                        actin.boundary_force_update(domainSize);
                    }
                } else {
                    // indicate that some of cad is not out of the boundary
                    // we need the flip method to move them into boundary
                    // need the smallest x and y position
                    double smallX = crossCluster.get(0).getBead_position()[0];
                    double smallY = crossCluster.get(0).getBead_position()[1];
                    for (Cadherin cad :crossCluster) {
                        if (Math.abs(cad.getBead_position()[0]) < Math.abs(smallX)) {
                            smallX = cad.getBead_position()[0];
                        }

                        if (Math.abs(cad.getBead_position()[1]) < Math.abs(smallY)) {
                            smallY = cad.getBead_position()[1];
                        }
                    }
                    // then we flip the corresponding cad
                    for (Cadherin cad : crossCluster) {
                        flip(cad, smallX, smallY);
                    }

                    // finally flip corresponing actin
                    for (Actin actin : actinCluster) {
                        flip(actin, smallX, smallY);
                    }
                }
            }
        }
    }

    private void findUpdate(List<Cadherin> mainsingles, List<Cadherin> helpsingles, ArrayList<HashSet<Integer>> clusters1, ArrayList<HashSet<Integer>> clusters2, int label, Map<Integer, Integer> info) {

        // this method (1) find the two-layer cluster (2) update the position
        double domainSize = 1;
        while (mainsingles.size() > 0) {// to ensure we always get the last cad in the singles
            boolean level = label != 1;

            // step 1: record the cluster + trans in two layers starting from 3.0 X, S in layer1
            ArrayList<Cadherin> crossCluster = new ArrayList<>(); // store the cadherin in this two-layer cluster
            ArrayList<Actin> actinCluster = new ArrayList<>(); // store the actins bind
            ArrayList<Integer> nextLevelCadindex = new ArrayList<>();
            ArrayList<Integer> curreentLevelCadindex = new ArrayList<>();
            boolean jump = true;
            int crossSize = 0;
            int i = mainsingles.size() - 1;
            Cadherin cad1 = mainsingles.get(i);
            crossCluster.add(cad1);
            //if (cad1.getActState() == 1) actinCluster.add(actinList_1.get(cad1.getActBoundIndex()));
            if (cad1.getActState() > 0) {
                if (cad1.getBead_position()[2] == 0) {
                    // cad1 is from layer1
                    actinCluster.add(actinList_1.get(cad1.getActBoundIndex()));
                } else {
                    // cad1 is from layer2
                    actinCluster.add(actinList_2.get(cad1.getActBoundIndex()));
                }
            }
            if (testBoundary(cad1, domainSize)) crossSize++;

            curreentLevelCadindex.add(cad1.getCadBoundIndex());
            mainsingles.remove(cad1); // delete current cad from the singles

            // jump between the first and second layer
            while (jump) {
                nextLevelCadindex.clear();
                for (int index : curreentLevelCadindex) {
                    if (level) {
                        // here we search the first clusters
                        Cadherin trBound = cadList_1.get(index);
                        if (!crossCluster.contains(trBound)) {
                            crossCluster.add(trBound);
                            if (trBound.getActState() > 0) actinCluster.add(actinList_1.get(trBound.getActBoundIndex()));
                            if (testBoundary(trBound, domainSize)) crossSize++;
                            for (int j = clusters1.size() - 1; j >= 0; j--) {
                                if (clusters1.get(j).contains(index)) {
                                    for (int cisindex : clusters1.get(j)) {
                                        if (cisindex != index) {
                                            Cadherin ciscad = cadList_1.get(cisindex);
                                            crossCluster.add(ciscad);
                                            if (ciscad.getActState() > 0) actinCluster.add(actinList_1.get(ciscad.getActBoundIndex()));
                                            if (testBoundary(ciscad, domainSize)) crossSize++;

                                            if (!ciscad.isCisM()) {
                                                int indexcad2 = ciscad.getCadBoundIndex();
                                                Cadherin cad2 = cadList_2.get(indexcad2);
                                                if (!crossCluster.contains(cad2)) {
                                                    if (cad2.getCisState() == 0) {
                                                        crossCluster.add(cad2);
                                                        if (cad2.getActState() > 0) actinCluster.add(actinList_2.get(cad2.getActBoundIndex()));
                                                        if (testBoundary(cad2, domainSize)) crossSize++;
                                                        // remove this cad2 from singles otherwise we repeat the work
                                                        if (label == 1) {
                                                            helpsingles.remove(cad2);
                                                        } else {
                                                            mainsingles.remove(cad2);
                                                        }
                                                    } else {
                                                        nextLevelCadindex.add(indexcad2);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    // remove this cluster in the cluster2
                                    clusters1.remove(clusters1.get(j));
                                    break;
                                }
                            }
                        }
                    } else {
                        // search the second clusters
                        Cadherin trBound = cadList_2.get(index);
                        if (!crossCluster.contains(trBound)) {

                            crossCluster.add(trBound);
                            if (trBound.getActState() > 0) actinCluster.add(actinList_2.get(trBound.getActBoundIndex()));
                            if (testBoundary(trBound, domainSize)) crossSize++;
                            for (int j = clusters2.size() - 1; j >= 0; j--) {

                                if (clusters2.get(j).contains(index)) {
                                    for (int cisindex : clusters2.get(j)) {
                                        if (cisindex != index) {
                                            Cadherin ciscad = cadList_2.get(cisindex);
                                            crossCluster.add(ciscad);
                                            if (ciscad.getActState() > 0) actinCluster.add(actinList_2.get(ciscad.getActBoundIndex()));
                                            if (testBoundary(ciscad, domainSize)) crossSize++;

                                            if (!ciscad.isCisM()) {
                                                int indexcad2 = ciscad.getCadBoundIndex();
                                                Cadherin cad2 = cadList_1.get(indexcad2);
                                                if (!crossCluster.contains(cad2)) {
                                                    if (cad2.getCisState() == 0) {
                                                        crossCluster.add(cad2);
                                                        if (cad2.getActState() > 0) actinCluster.add(actinList_1.get(cad2.getActBoundIndex()));
                                                        if (testBoundary(cad2, domainSize)) crossSize++;
                                                        // remove this cad2 from singles otherwise we repeat the work
                                                        if (label == 1) {
                                                            mainsingles.remove(cad2);
                                                        } else {
                                                            helpsingles.remove(cad2);
                                                        }
                                                    } else {
                                                        nextLevelCadindex.add(indexcad2);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    // remove this cluster in the cluster2
                                    clusters2.remove(clusters2.get(j));
                                    break;
                                }
                            }
                        }
                    }
                }
                // here we iterate through the cluster contains the index in the currentlevel
                // we update several condtions
                if (nextLevelCadindex.size() == 0) {
                    jump = false;
                } else {
                    curreentLevelCadindex = (ArrayList<Integer>)(nextLevelCadindex.clone());
                    level = !level; // jump to the next level
                }
            }

            // count the corresponding two layer cluster size
            int twolayer_size = crossCluster.size();
            if (!info.containsKey(twolayer_size)) {
                info.put(twolayer_size, 1);
            } else {
                info.put(twolayer_size, info.get(twolayer_size) + 1);
            }

            // here we find the corresponding two level cluster
            // decide flip or not
            // we only move the cad and actin
            if (crossSize > 0) {
                if (crossSize == crossCluster.size()) {
                    // this indicates all cad in this cluster out of the boundary
                    // under this situation, we just do the normal boundary update
                    for (Cadherin cad : crossCluster) {
                        cad.boundary_force_update(domainSize);
                    }
                    for (Actin actin : actinCluster) {
                        actin.boundary_force_update(domainSize);
                    }
                } else {
                    // indicate that some of cad is not out of the boundary
                    // we need the flip method to move them into boundary
                    // need the smallest x and y position
                    double smallX = crossCluster.get(0).getBead_position()[0];
                    double smallY = crossCluster.get(0).getBead_position()[1];
                    for (Cadherin cad :crossCluster) {
                        if (Math.abs(cad.getBead_position()[0]) < Math.abs(smallX)) {
                            smallX = cad.getBead_position()[0];
                        }

                        if (Math.abs(cad.getBead_position()[1]) < Math.abs(smallY)) {
                            smallY = cad.getBead_position()[1];
                        }
                    }
                    // then we flip the corresponding cad
                    for (Cadherin cad : crossCluster) {
                        flip(cad, smallX, smallY);
                    }

                    // finally flip corresponing actin
                    for (Actin actin : actinCluster) {
                        flip(actin, smallX, smallY);
                    }
                }
            }
        }
    }

    private void flip(Actin actin, double smallX, double smallY) {

        // flip the actin based on the smallX and smallY, we just want to make sure the relative distance
        if (smallX > 0) {
            actin.setBead_position_x(smallX - Math.abs(actin.getBead_position()[0] - smallX));
        } else {
            actin.setBead_position_x(smallX + Math.abs(actin.getBead_position()[0] - smallX));
        }

        // update the y position
        if (smallY > 0) {
            actin.setBead_position_y(smallY - Math.abs(actin.getBead_position()[1] - smallY));
        } else {
            actin.setBead_position_y(smallY + Math.abs(actin.getBead_position()[1] - smallY));
        }
    }

    private void flip(Cadherin cad, double smallX, double smallY) {

        // this method is to flip the bead position of cad to make sure it is in the boundary
        // update the x position
        if (smallX > 0) {
            cad.setBead_position_x(smallX - Math.abs(cad.getBead_position()[0] - smallX));
        } else {
            cad.setBead_position_x(smallX + Math.abs(cad.getBead_position()[0] - smallX));
        }

        // update the y position
        if (smallY > 0) {
            cad.setBead_position_y(smallY - Math.abs(cad.getBead_position()[1] - smallY));
        } else {
            cad.setBead_position_y(smallY + Math.abs(cad.getBead_position()[1] - smallY));
        }
    }

    private boolean testBoundary(Cadherin cad, double domainSize) {
        return (cad.getBead_position()[0] > domainSize/2 || cad.getBead_position()[0] < -domainSize/2 ||
                cad.getBead_position()[1] > domainSize/2 || cad.getBead_position()[1] < -domainSize/2);
    }

    private ArrayList<HashSet<Integer>> moveCisCad2(ArrayList<Cadherin> cadList, ArrayList<Actin> actinList, double force) {

        double actinFrictional = 0.108;
        double actink_it = 2e3;
        // in this method we move the cad in cis according to the small cluster
        // the container to store all small cluster in the corresponding layer
        ArrayList<HashSet<Integer>> Clusters = new ArrayList<>(); // only for cases cisstate = 1 or 2
        ArrayList<Integer> clusterindex = new ArrayList<>();
        Clusters.add(new HashSet<>()); // the first one is empty.
        // directly filter out all singles with either weakly or strongly bond
        List<Cadherin> Singles = cadList.stream().filter(cad -> (cad.getState() == 3 && cad.getCisState() == 0 && cad.getActState() > 0)).collect(Collectors.toList());
        // then we can focus on the cisState = 1
        // after we filter out the clusters with the head 3.1, the rest will be the head 3.2
        for (Cadherin cad : cadList) {
            if (cad.getState() == 3 && cad.getCisState() == 1) {
                int index = cadList.indexOf(cad);
                if (!clusterindex.contains(index)) {
                    HashSet<Integer> newCluster = new HashSet<>();
                    // all previous clusters do not have this cad, therefore we can use this cad as the head to start a new cluster
                    // when the head of cluster is in the state of 3.1
                    findCluster1(newCluster, cad, cadList);
                    Clusters.add(newCluster);
                    clusterindex.addAll(newCluster);
                }
            }
        }

        // we deal with the rest cluster starting with head 3.2
        for (Cadherin cad : cadList) {
            if (cad.getState() == 3 && cad.getCisState() == 2) {
                int index = cadList.indexOf(cad);
                if (!clusterindex.contains(index)) {
                    HashSet<Integer> newCluster = new HashSet<>();
                    findCluster2(newCluster, cad, cadList);
                    Clusters.add(newCluster);
                    clusterindex.addAll(newCluster);
                }
            }
        }

        // then we can move them accordingly
        // double force = 0;
        // move the cad (3,0)
        for (Cadherin cad : Singles) {
            double theta = -Math.PI + Math.random()*2*Math.PI; // generate a random direction
            Actin actin = actinList.get(cad.getActBoundIndex());
            actinMove(actin, force, theta, actinFrictional); // move the actin
            double caforce = distance(cad, actin) * actink_it;
            // determine the break or not according to this force
            // get the state of actin
            int actinState = cad.getActState();
            double kub = deterKub(cad, actinState, caforce);
            if (Math.random() <= kub*timestep) {
                // this indicates we need to remove the actin
                cad.setActState(0);
                actinList.set(cad.getActBoundIndex(), null);
            } else {
                boundForceUpdate(cad, actin);
                actin.setForce(actink_it * distance(cad, actin));
            }
        }

        // move the (3.1) and (3.2) cluster
        for (HashSet<Integer> cluster : Clusters) {
            int boundactsize = 0; // cadherin bound with actin
            double totalSpforce = 0;
            double theta = -Math.PI + Math.random()*2*Math.PI;
            // count how many cadherin in the cluster has actin bound
            for (int index : cluster) {
                if (cadList.get(index).getActState() > 0) boundactsize++;
            }
            // distribute the force and move
            if (boundactsize > 0) { // we have the actin binding for this clsuter
                // 1. distribute the force
                double disforce = force / boundactsize;
                for (int index : cluster) {
                    // get corresponding cad
                    Cadherin cad = cadList.get(index); // get corresponding cadherin
                    int actState = cad.getActState();
                    if (actState > 0) { // we get the corresponding binding actin
                        // means we have the corresponding bond
                        Actin actin = actinList.get(cad.getActBoundIndex());
                        // move the actin with the distribute force
                        actinMove(actin, disforce, theta, actinFrictional);
                        double caforce = actink_it * distance(cad, actin);
                        double kub = deterKub(cad, actState, caforce);
                        // determine the break or not
                        if (Math.random() <= kub * timestep) { // break the bond
                            cad.setActState(0);
                            actinList.set(cad.getActBoundIndex(), null);
                        } else {
                            // the bond is bot break
                            totalSpforce += caforce; // this is the spring force beween actin and cadherin
                        }
                    }
                }
            }

            // distribute the total spring force and move each cadherin in the cluster
            double disSpforce = totalSpforce / cluster.size(); // even if we have an empty hashset, it doesn't matter
            // step 3: move the cad
            for (int index : cluster) {
                Cadherin cad = cadList.get(index);
                cadMove(cad, disSpforce, theta);
                if (cad.getActState() > 0) {
                    Actin actin = actinList.get(cad.getActBoundIndex());
                    actin.setForce(actink_it * distance(cad, actin));
                }
            }
        }
        return Clusters;
    }

    private double deterKub(Cadherin cad, int actinState, double caforce) {

        // to decide the kub of bond between cadherin and actin
        if (actinState == 1) {
            // determine whether switch to the strongly bond or not
            if (Math.random() <= f12(caforce) * timestep) {
                cad.setActState(2); // change the state to 2
                System.out.println("Switch the bond state from weak state 1 to strong state 2!");
                return fstrong(caforce);
            } else {
                return fweak(caforce);
            }
        } else {
            // determine whether switch to the weakly bound and determine the kub accordingly
            if (Math.random() <= f21(caforce) * timestep) {
                cad.setActState(1);
                System.out.println("Switch the bond state from strong state 2 to weak state 1!");
                return fweak(caforce);
            } else {
                return fstrong(caforce);
            }
        }
    }

    private void findCluster2(HashSet<Integer> cluster, Cadherin head, ArrayList<Cadherin> cadList) {

        // here we find the cluster starting with cad in the state2
        // if there is any possibility that the 3.2 is before 3.1 ??
        int headindex = cadList.indexOf(head);
        cluster.add(headindex);
        int previndex = headindex;
        int currentindex = head.getCadCisBoundIndex().get(0);
        Cadherin current = cadList.get(currentindex);
        int endindex = head.getCadCisBoundIndex().get(1);
        int nextindex = -1;

        // use the while loop to find the cluster
        while (nextindex != endindex) {
            // find the next index
            for (int index : current.getCadCisBoundIndex()) {
                if (index != previndex) {
                    nextindex = index;
                }
            }
            cluster.add(currentindex);
            previndex = currentindex;
            currentindex = nextindex;
            try {
                current = cadList.get(currentindex);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        cluster.add(endindex);
    }

    private void findCluster1(HashSet<Integer> cluster, Cadherin head, ArrayList<Cadherin> cadList) {

        int headindex = cadList.indexOf(head);
        cluster.add(headindex); // indicating this cad is the head of the cluster
        int previndex = headindex;
        int currentindex = 0; // check
        try {
            currentindex = head.getCadCisBoundIndex().get(0);
        } catch (Exception e) {
            e.printStackTrace();
        }
        Cadherin current = cadList.get(currentindex);
        int nextindex = -1;

        // use the while loop
        while (current.getCisState() == 2) {
            // find the nextindex
            for (int integer : current.getCadCisBoundIndex()) {
                if (integer != previndex) {
                    nextindex = integer;
                }
            }
            // add the current cad in the cluster
            cluster.add(currentindex);
            // change the index
            previndex = currentindex;
            currentindex = nextindex;
            current = cadList.get(currentindex);
        }

        // solve the tail of cluster
        if (current.getCisState() == 1 && current.getCadCisBoundIndex().get(0) == previndex) {
            cluster.add(currentindex);
        }
    }

    private void cadMove(Cadherin cad, double distributedForce, double theta) {
        double[] bead_force = {distributedForce*Math.cos(theta), distributedForce*Math.sin(theta)};
        for (int i = 0; i < 2; i++) {
            cad.getBead_position()[i] -= bead_force[i]*timestep/cad.getFrictional();
        }
    }

    private void deterFate(Cadherin cad, Cadherin cad2, Actin actin1, Actin actin2, double actink_it, double k_it) {

        double cutoff_2 = 0.001; // unit: micrometers
        double cutoff = 0.003;
        double k_fxs = 86; // x -> s unit: 1/s
        double k_bxs = 0.8;
        double p_fxs = k_fxs * timestep;
        double p_bxs = k_bxs * timestep;
        double dist = distance(cad, cad2); // sometimes, this dist is zero, therefore thinking about the boundary condition
        double force = k_it*dist;
        // decide it is in X or S dimer in the cis
        if (cad.isCisX() /*&& !(cad.isCisM())*/) { // in the X-dimer cad.isCisM() can be deleted
            double kub = xCatch(force);
            if (Math.random() <= kub*timestep) {

                cad.setCisM(true); cad.setCisX(false); cad.setCadBoundIndex(-1); cad.setForce(0);
                cad2.setCisM(true); cad2.setCisX(false); cad2.setCadBoundIndex(-1); cad2.setForce(0);
                // new boundary condition!
                if (cad.getState() == 3 && cad.getCisState() == 0) {
                    cad.setState(0);
                    if (cad.getActState() > 0) {
                        actinList_1.set(cad.getActBoundIndex(), null);
                        cad.setActState(0);
                    }
                    cad.setCisM(false);
                }
                if (cad2.getState() == 3 && cad2.getCisState() == 0) {
                    cad2.setState(0);
                    if (cad2.getActState() > 0) {
                        actinList_2.set(cad2.getActBoundIndex(), null);
                        cad2.setActState(0);
                    }
                    cad2.setCisM(false);
                }
            } else {
                // the bond does not break
                boundForceUpdate(cad, cad2/*, dist*/); // update the position
                double ndist = distance(cad, cad2);
                double nforce = k_it * ndist;
                // update the force if they are in X or S
                cad.setForce(nforce);
                cad2.setForce(nforce);

                if (ndist <= cutoff_2 && Math.random() <= p_fxs) {
                    cad.setCisX(false);
                    cad2.setCisX(false);
                    // update the corresponding actin
                }

                if (null != actin1) {
                    boundForceUpdate(cad, actin1/*, distance(cad, actin1)*/);
                    actin1.setForce(actink_it * distance(cad, actin1));
                }
                if (null != actin2) {
                    boundForceUpdate(cad2, actin2/*, distance(cad2, actin2)*/);
                    actin2.setForce(actink_it * distance(cad2, actin2));
                }
                return;
            }
        } else { // still we don't need the cisM condition
            // as S-dimer in the cis
            double kub = sSlip(force);
            if (Math.random() <= kub*timestep) {
                cad.setCisM(true); cad.setCisX(false); cad.setCadBoundIndex(-1); cad.setForce(0); // setCisX redundant
                cad2.setCisM(true); cad2.setCisX(false); cad2.setCadBoundIndex(-1); cad2.setForce(0);
                if (cad.getState() == 3 && cad.getCisState() == 0) {
                    cad.setState(0);
                    if (cad.getActState() > 0) {
                        actinList_1.set(cad.getActBoundIndex(), null);
                        cad.setActState(0);
                    }
                    cad.setCisM(false);
                }
                if (cad2.getState() == 3 && cad2.getCisState() == 0) {
                    cad2.setState(0);
                    if (cad2.getActState() > 0) {
                        actinList_2.set(cad2.getActBoundIndex(), null);
                        cad2.setActState(0);
                    }
                    cad2.setCisM(false);
                }
            } else {
                boundForceUpdate(cad, cad2/*, dist*/);
                double ndist = distance(cad, cad2);
                double nforce = k_it*ndist;
                cad.setForce(nforce);
                cad2.setForce(nforce);


                // update the corresponding actin
                if (null != actin1) {
                    boundForceUpdate(cad, actin1);
                    actin1.setForce(actink_it * distance(cad, actin1));
                }
                if (null != actin2) {
                    boundForceUpdate(cad2, actin2);
                    actin2.setForce(actink_it * distance(cad2, actin2));
                }
                if (Math.random() <= p_bxs) {
                    // from S-dimer in the cis back to X dimer
                    cad.setCisX(true);
                    cad2.setCisX(true);
                }
            }
        }
    }

    private void actinMove(Actin actin1, double force1, double theta1, double actinFrictional) {

        // move the actin in the corresoonding x and y direction
        double[] bead_force = {force1*Math.cos(theta1), force1*Math.sin(theta1)};
        for (int i = 0; i < 2; i++) {
            actin1.getBead_position()[i] -= bead_force[i]*timestep/actinFrictional;
        }
    }
}