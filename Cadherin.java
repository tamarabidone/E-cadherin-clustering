package com.company;

import java.util.*;
import java.lang.Math; // to get the constant value pi

/**
 * This module is the first version cadherin object with corresponding features
 * Here we only focus on the thermal and boundary force, but leave the actin cortical force in the other module
 * Basically, we describe the feature of cadherin in the space and time.
 * @author yangchen
 */

public class Cadherin {
    // instance variables
    private double time = 0.0;
    private double [] bead_position;
    private double [] ref_position; // todo: delete no need for this line
    private double [] bead_force; // here I assume the force only change at x, y plane
    private double frictional; // frictional coefficient, D = (1/frictional) * kT

    // record the corresponding state of cadherin
    private int state = 0; // to indicate the state of the cadherin, can have three different states, 0-M, 1-X, 2-S, 3-cis
    private int cisState = 0; // to indicate the number of cis-interface formed
    private int actState = 0; // if we allow the actState = 1 -> weakly bound, actState = 2 means strongly bound

    // to record the force between the dimers
    private double force = 0;

    // record the corresponding bounded trans cad or cis cad or actin
    private int cadBoundIndex = -1; // the index of bound Cadherin
    private ArrayList<Integer> cadCisBoundIndex = new ArrayList<Integer>(2); // to store the list of index of cis-interface bound Cadherin
    private int actBoundIndex = -1;

    // indicate whether it is cisM or cisX
    // cooperate with the State = 3 to indicate the state
    private boolean cisM = false;
    private boolean cisX = false;

    // constructor
    public Cadherin(double [] ini_position, double frictional) {
        this.ref_position = ini_position.clone(); // todo: delete this line no need for this project
        this.bead_position = ini_position.clone();
        this.frictional = frictional;
    }

    // The thermal force and the boundary force part are the stochastic force to update the position of the bead
    // which is the first part of the stochastic diffusion.

    public void thermal_force_update(double temperature, double dt) {
        bead_force = new double[2]; // each time renew the bead_force
        //this.time += dt;
        Random rand = new Random();
        double thermal_force = Math.sqrt(2*0.0000138*temperature*frictional/dt);

        // for all state this step should be the same
        // they all need to follow the thermal force
        // to calculate the thermal force and update the bead_force
        bead_force[0] += thermal_force*rand.nextGaussian();
        bead_force[1] += thermal_force*rand.nextGaussian();

        move(dt);

        // here I don't consider the thermal force in the z direction yet.
    }

    public void boundary_force_update(double domainSize) {
        // update the x direction
        if (bead_position[0] > domainSize/2) {
            bead_position[0] -= domainSize;
        }
        if (bead_position[0] < -domainSize/2) {
            bead_position[0] += domainSize;
        }
        // update the y direction
        if (bead_position[1] > domainSize/2) {
            bead_position[1] -= domainSize;
        }
        if (bead_position[1] < -domainSize/2) {
            bead_position[1] += domainSize;
        }
    }

    // can make this into public method
    private void move(double dt) {
        for (int i = 0; i < 2; i++) {
            bead_position[i] -= bead_force[i]*dt/frictional;
        }
    }

    @Override
    public String toString() {
        return "Cadherin{" +
                "time=" + time +
                ", bead_position=" + Arrays.toString(bead_position) +
                ", bead_force=" + Arrays.toString(bead_force) +
                ", state=" + state +
                ", cisState=" + cisState +
                ", actState=" + actState +
                ", force=" + force +
                ", cadBoundIndex=" + cadBoundIndex +
                ", cadCisBoundIndex=" + cadCisBoundIndex +
                ", actBoundIndex=" + actBoundIndex +
                ", cisM=" + cisM +
                ", cisX=" + cisX +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Cadherin cadherin = (Cadherin) o;
        return Double.compare(cadherin.time, time) == 0 &&
                state == cadherin.state &&
                cisState == cadherin.cisState &&
                actState == cadherin.actState &&
                Double.compare(cadherin.force, force) == 0 &&
                cadBoundIndex == cadherin.cadBoundIndex &&
                actBoundIndex == cadherin.actBoundIndex &&
                cisM == cadherin.cisM &&
                cisX == cadherin.cisX &&
                Arrays.equals(bead_position, cadherin.bead_position) &&
                cadCisBoundIndex.equals(cadherin.cadCisBoundIndex);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(time, state, cisState, actState, force, cadBoundIndex, cadCisBoundIndex, actBoundIndex, cisM, cisX);
        result = 31 * result + Arrays.hashCode(bead_position);
        return result;
    }

    // set the value of the cadBoundIndex
    public void setCadBoundIndex(int value) { this.cadBoundIndex = value; }

    // set the value of the state of the Cadherin
    public void setState(int value) { this.state = value; }

    public void setForce(double value) { this.force = value; }

    // set the bead_position in x axis
    public void setBead_position_x(double position_x) { this.bead_position[0] = position_x; }

    public void setBead_position_y(double position_y) { this.bead_position[1] = position_y; }

    // method the set the time
    public void setTime(double dt) { this.time += dt; }

    // method to add index of cisbound cad
    public void setCadCisBoundIndex(int index) { cadCisBoundIndex.add(index); }

    // set the value of cisState
    public void setCisState(int state) { this.cisState = state; }

    public void setActState(int state) { this.actState = state; }


    public void setCisM(boolean cisM) {
        this.cisM = cisM;
    }

    public void setCisX(boolean cisX) {
        this.cisX = cisX;
    }

    public void setActBoundIndex(int actBoundIndex) {
        this.actBoundIndex = actBoundIndex;
    }

    // set the bead_force
    //  public void setBead_force(double[] bead_force) { this.bead_force = bead_force.clone(); }

    // get the bead_position for the user
    public double[] getBead_position() {
        return this.bead_position;
    }

    // the following two getters are for the move together when the state == 1 / 2
    // get the bead_force for the user
    public double[] getBead_force() { return this.bead_force; }

    // getter of time for the user
    public double getTime() {
        return this.time;
    }

    public double[] getRef_position() {return this.ref_position; }

    public double getFrictional() {
        return frictional;
    }

    // get the current state of cadherin
    public int getState() { return this.state; }

    // get the index of bounded cadherin
    public int getCadBoundIndex() { return this.cadBoundIndex; }

    // get the bound force between the cadherin
    public double getForce() { return this.force; }

    // get the bound index of cis-interface cadherin
    public ArrayList<Integer> getCadCisBoundIndex() { return this.cadCisBoundIndex; }

    // get the cisState of cad
    public int getCisState() { return this.cisState; }

    // get the actState
    public int getActState() { return this.actState; }

    //public boolean isIsadd() { return isadd; }

    public int getActBoundIndex() {
        return actBoundIndex;
    }

    public boolean isCisM() {
        return cisM;
    }

    public boolean isCisX() {
        return cisX;
    }
}
