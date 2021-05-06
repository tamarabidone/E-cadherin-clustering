package com.company;

import java.util.Arrays;

/**    This class file is model the F-actin filament in the project.
 *     The following are the assumptions we think for the F-actin
 *     1. We still model the F-actin as a bead on another surface to form the catch bound with the Cadherin in the
 *        corresponding layer
 *     2. The binding between the F-actin is one-to-one once the state of cadherin is 3 (forward rate constant?) MAYBE PEOPLE CAN TEST THE CONCENTRATION OF ACTIN
 *     3. The unbinding will happen when we apply certain force on the F-actin bead to move their bead-position on purpose
 *     4. We also assume that the distance between cadherin and F-actin (Todo: find ref for this distance?)
 * @author yang chen
 */

public class Actin {
    // Some variable we wanna store at this stage
    private double[] bead_position; // same as the cadherin position except the z-position
    private int BoundcadIndex = -1;
    private double frictional = 0.108; // unit pN * sec / um
    private double force; // use later to break the bound. Todo: combine with the catch bond type (rate constant?)

    public Actin(double[] bead_position) {
        this.bead_position = bead_position.clone();
    }

    // setter the change the binding index
    public void setBoundcadIndex(int index) { this.BoundcadIndex = index; }


    public void setForce(double force) {
        this.force = force;
    }

    public double[] getBead_position() {
        return bead_position;
    }

    public double getFrictional() {
        return frictional;
    }

    public double getForce() {
        return force;
    }

    public void setBead_position_x(double v) {
        this.bead_position[0] = v;
    }

    public void setBead_position_y(double v) {
        this.bead_position[1] = v;
    }

    @Override
    public String toString() {
        return "Actin{" +
                "bead_position=" + Arrays.toString(bead_position) +
                ", BoundcadIndex=" + BoundcadIndex +
                ", force=" + force +
                '}';
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
}
