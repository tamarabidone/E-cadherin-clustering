package com.company;

import java.io.FileWriter; // to write the results to the file
import java.io.IOException;
import java.util.ArrayList;
import java.lang.Math;
import java.util.Arrays;
import java.util.Map;

/**
 * This is the main module to test the Cadherin object and simulate the force interaction between
 * two layers of Cadherin
 * @author yangchen
 */

public class Main {

    // function to calculate the mean-squared displacement
    public static double msd(double[] ref_pos, double[] curr_pos) {
        return Math.pow(curr_pos[0] - ref_pos[0], 2) + Math.pow(curr_pos[1] - ref_pos[1],2);
    }

    // function to test whether the cad is out of the domain or not
    public static boolean testBoundary(Cadherin cad, double domainSize) {
        return (cad.getBead_position()[0] > domainSize/2 || cad.getBead_position()[0] < -domainSize/2 ||
                cad.getBead_position()[1] > domainSize/2 || cad.getBead_position()[1] < -domainSize/2);
    }

    // function to Relative Boundary Condition update
    public static void relativeBFU(Cadherin cad1, Cadherin cad2, double domainSize) {
        // here treat cad1 as the ref cad
        // ref cad means that this cad is outside the boundary, therefore update the position of cad2 relatively
        if (cad1.getBead_position()[0] > domainSize/2) {
            double update_pos_x = cad1.getBead_position()[0]+ Math.abs(cad1.getBead_position()[0]-cad2.getBead_position()[0]);
            cad2.setBead_position_x(update_pos_x);
        }

        if (cad1.getBead_position()[0] < -domainSize/2) {
            double update_pos_x = cad1.getBead_position()[0]- Math.abs(cad1.getBead_position()[0]-cad2.getBead_position()[0]);
            cad2.setBead_position_x(update_pos_x);
        }

        // y direction
        if (cad1.getBead_position()[1] > domainSize/2) {
            double update_pos_y = cad1.getBead_position()[1]+ Math.abs(cad1.getBead_position()[1]-cad2.getBead_position()[1]);
            cad2.setBead_position_y(update_pos_y);
        }

        if (cad1.getBead_position()[1] < -domainSize/2) {
            double update_pos_y = cad1.getBead_position()[1]- Math.abs(cad1.getBead_position()[1]-cad2.getBead_position()[1]);
            cad2.setBead_position_y(update_pos_y);
        }
    }


    // function to relatively move the thermal of the cadherin in the state = 2
    public static void relativeTFU(Cadherin cad1, Cadherin cad2, double frictional, double dt) {
        // still the cad1 is the ref cad
        // when call this function the ref cad already got updated thermally
        double[] frc = cad1.getBead_force().clone();
        double update_pos_x, update_pos_y;
        update_pos_x = cad2.getBead_position()[0] - frc[0]*dt/frictional;
        cad2.setBead_position_x(update_pos_x);
        update_pos_y = cad2.getBead_position()[1] - frc[1]*dt/frictional;
        cad2.setBead_position_y(update_pos_y);
    }


    // function to update boundary accordingly
    public static void update_Boundary(Cadherin cad1, Cadherin cad2, double domainSize) {
        if (testBoundary(cad1, domainSize)) {
            // here already update the ref cad position
            cad1.boundary_force_update(domainSize);
            if (testBoundary(cad2, domainSize)) {
                cad2.boundary_force_update(domainSize);
            } else {// cad2 is inside the domain
                relativeBFU(cad1, cad2, domainSize);
            }
        } else { // cad1 is inside the domain
            if (testBoundary(cad2, domainSize)) {
                cad2.boundary_force_update(domainSize);
                relativeBFU(cad2, cad1, domainSize);
            }
        }
    }

    // to calculate the proportion of at each record time point
    public static double[] count(ArrayList<Cadherin> cadList, String filename) {
        int[] counts = new int[4]; //  index0 = M , index1 = X, index2 = S, index3 = cis
        int[] total = new int[3]; // all M, X, S in the status of free or cis cluster
        int[] free = new int[3]; // all M, X, S not in the cis
        int[] incis = new int[3]; // all M, X, S in the cis
        int[] forcefree = new int[2]; // only X and S experiencing force 0 X 1 S
        int[] forceincis = new int[2]; // still only X and S experiencign the fore
        String path = "./" + filename + "_prop.csv";
        String path2 = "./" + filename + "_free"  + ".csv";
        String path3 = "./" + filename + "_incis" + ".csv";
        String path4 = "./" + filename + "_total" + ".csv";


        for (Cadherin cad : cadList) { // to record all states
            if (cad.getState() == 3) {
                counts[3]++; // this is in the cis cluster
                if (cad.isCisM()) {
                    total[0]++; // the total M ++
                    incis[0]++; // this M is in the cis
                } else if (cad.isCisX()) {
                    forceincis[0] += cad.getForce(); // force for X in the cis
                    total[1]++; // the total X ++
                    incis[1]++; // the X is in the cis
                } else {
                    forceincis[1] += cad.getForce(); // force for S in the cis
                    total[2]++; // the total S ++
                    incis[2]++; // the S is in the cis
                }
            } else if (cad.getState() == 2) {
                forcefree[1] += cad.getForce();
                counts[2]++; // the S dimer
                total[2]++; // the S dimer total
                free[2]++; // the free S dimer
            } else if (cad.getState() == 1) {
                forcefree[0] += cad.getForce();
                counts[1]++; // the X dimer
                total[1]++;
                free[1]++;
            } else {
                counts[0]++;
                total[0]++;
                free[0]++;
            }
        }

        double[] propotion = new double[4];
        for (int i= 0; i < propotion.length; i++) {
            propotion[i] = counts[i]*1.0 / cadList.size();
        }

        double[] totalprop = new double[3];
        double[] freeprop = new double[3];
        double[] cisprop = new double[3];
        for (int i = 0; i < totalprop.length; i++) {
            totalprop[i] = total[i]*1.0 / cadList.size();
            freeprop[i] = free[i]*1.0 / cadList.size();
            cisprop[i] = incis[i]*1.0 / cadList.size();
        }

        double[] avgforcecis = new double[2];
        double[] avgforcefree = new double[2];
        // calculate the average force
        for (int i = 0; i < avgforcecis.length; i++) {
            if (incis[i+1] != 0) {
                avgforcecis[i] = forceincis[i]*1.0 / incis[i+1];
            } else {
                avgforcecis[i] = 0;
            }
            if (free[i+1] != 0) {
                avgforcefree[i] = forcefree[i]*1.0 / free[i+1];
            } else {
                avgforcefree[i] = 0;
            }
        }

        FileWriter fw = null; // write the normal prop
        FileWriter fw2 = null; // write the free prop + free force
        FileWriter fw3 = null; // write the cis prop + cisave force
        FileWriter fw4 = null; // write the total prop
        try {
            fw = new FileWriter(path, true);
            fw2 = new FileWriter(path2, true);
            fw3 = new FileWriter(path3, true);
            fw4 = new FileWriter(path4, true);
            fw.write(String.format("%6.4f," + "%6.4f," + "%6.4f," + "%6.4f" + "\n", propotion[0], propotion[1], propotion[2], propotion[3]));
            fw.flush();
            fw2.write(String.format("%6.4f," + "%6.4f," + "%6.4f," + "%6.4f," + "%6.4f" + "\n", freeprop[0], freeprop[1], freeprop[2], avgforcefree[0], avgforcefree[1]));
            fw2.flush();
            fw3.write(String.format("%6.4f," + "%6.4f," + "%6.4f," + "%6.4f," + "%6.4f" + "\n", cisprop[0], cisprop[1], cisprop[2], avgforcecis[0], avgforcecis[1]));
            fw3.flush();
            fw4.write(String.format("%6.4f," + "%6.4f," + "%6.4f" + "\n", totalprop[0], totalprop[1], totalprop[2]));
            fw4.flush();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (null != fw) {
                try {
                    fw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            if (null != fw2) {
                try {
                    fw2.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            if (null != fw3) {
                try {
                    fw3.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            if (null != fw4) {
                try {
                    fw4.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        return propotion;
    }


    // function to simulate the binding between two layers
    public static void simulate(ArrayList<Cadherin> cadList1, ArrayList<Cadherin> cadList2, double totalTime,
                                Force_Interaction force_interaction, double timestep, String nameofFile,
                                double temperature, double domainSize, double frictional) {
        int size1 = cadList1.size(); int size2 = cadList2.size();
        String path = "./" + nameofFile + ".csv";
        String path2 = "./" + nameofFile + "clusterinfo" + ".txt";
        String nameofFile2 = nameofFile + "_1s";
        String path3 = "./" + nameofFile2 + ".csv"; // to record what will happen in 1 second, with 0.1 step
        String path4 = "./" + nameofFile2 + "clsuterinfo" + ".txt";


        try {
            FileWriter fw = new FileWriter(path, true); // The full dataset file writer
            FileWriter fw5 = new FileWriter(path2, true); // store the info of clusters
            FileWriter fw6 = new FileWriter(path3, true); // 1s full dataset
            FileWriter fw7 = new FileWriter(path4, true); // 1s info of clusters

            fw.write(String.format("index," + "time," + "state," + "cisState," + "actState," + "iscisM," + "iscisX," + "X," + "Y," + "Z," + "BCI," + "FORCE," + "CisBCI1," + "CisBCI2" + "\n"));
            fw6.write(String.format("index," + "time," + "state," + "cisState," + "actState," + "iscisM," + "iscisX," + "X," + "Y," + "Z," + "BCI," + "FORCE," + "CisBCI1," + "CisBCI2" + "\n"));

            // simulate through
            for (int i = 0; i <= (int)(totalTime/timestep); i++) {
                //System.out.println("The " + i + "th interation!");
                for (Cadherin cad1 : cadList1) { cad1.setTime(timestep); }
                for (Cadherin cad2 : cadList2) { cad2.setTime(timestep); }
                // step 1: do the thermal and boundary update according to the state
                // First update the all Cadherin in two layers with state = 0
                // Second, use the Cadherin in layer 1 as the ref to update the thermal and boundary accordingly
                for (Cadherin cad1 : cadList1) {
                    if (cad1.getState() == 0) {
                        cad1.thermal_force_update(temperature, timestep);
                        cad1.boundary_force_update(domainSize);
                    }

                    if (cad1.getState() == 1) {
                        Cadherin BoundCad = cadList2.get(cad1.getCadBoundIndex());
                        // update thermal simultaneously (independently)
                        cad1.thermal_force_update(temperature, timestep);
                        BoundCad.thermal_force_update(temperature, timestep);

                        // update the boundary condition accordingly
                        update_Boundary(cad1, BoundCad, domainSize);
                    } // in this way all Cadherins in state = 1 get updated both thermally and boundary condition

                    if (cad1.getState() == 2) {
                        Cadherin BoundCad = cadList2.get(cad1.getCadBoundIndex());
                        // update the thermal of them
                        cad1.thermal_force_update(temperature, timestep);
                        relativeTFU(cad1, BoundCad, frictional, timestep);

                        // update the boundary condition accordingly
                        update_Boundary(cad1, BoundCad, domainSize);
                    }

                }

                for (Cadherin cad2 : cadList2) {
                    if (cad2.getState() == 0) {
                        cad2.thermal_force_update(temperature, timestep);
                        cad2.boundary_force_update(domainSize);
                    }
                }

                // step 2: force interaction update
                force_interaction.trans_interaction();

                // step 3: cis interaction update
                // Test the running time of this section
                long startTime = System.nanoTime();
                force_interaction.cis_interaction3();
                long endTime = System.nanoTime();
                long timeElapsed = endTime - startTime;
                if ((i % 100000) == 0) System.out.println("The running time of this cis_section is: " + timeElapsed);

                // step 3: actin interaction update
                ArrayList<Map<Integer, Integer>> info;
                info = force_interaction.actin_interaction3();

                // step 4: trans_interaction in the cis-cluster
                force_interaction.trans_cis_interaction();

                // step 5: record the state
                // record every 100000 iterations
                if ( (i% 100000)  == 0) {

                    System.out.println("-------------------------------------------------");
                    System.out.println("record!");
                    double[] count1 = count(cadList1, nameofFile);
                    //double[] count2 = count(cadList2, "actin200prosimf230_100_5");
                    //double[] count2 = count(cadList2);
                    System.out.println("the proportion for each state for layer one is: " + Arrays.toString(count1));
                    //System.out.println("the proportion for each state for layer tow is: " + Arrays.toString(count2));
                    long startTime1 = System.nanoTime();
                    for (Cadherin cad1 : cadList1) {
                        fw.write(String.format("layer1-%d," +  "%6.4f," + "%d," + "%d," + "%d," + "%b," + "%b," + "%6.4f," + "%6.4f," + "%6.4f," + "%d," + "%6.4f," + cad1.getCadCisBoundIndex().toString() + "," +"\n",
                                cadList1.indexOf(cad1),  cad1.getTime(), cad1.getState(), cad1.getCisState(), cad1.getActState(), cad1.isCisM(), cad1.isCisX(), cad1.getBead_position()[0], cad1.getBead_position()[1],
                                cad1.getBead_position()[2], cad1.getCadBoundIndex(), cad1.getForce()));
                    }
                    long endTime1 = System.nanoTime();
                    long timeElapsed1 = endTime1 - startTime1;
                    System.out.println("The running time of this file_writing is: " + timeElapsed1);

                    for (Cadherin cad2 : cadList2) {
                        fw.write(String.format("layer2-%d," + "%6.4f," + "%d," + "%d," + "%d," + "%b," + "%b,"+ "%6.4f," + "%6.4f," + "%6.4f," + "%d," + "%6.4f," + cad2.getCadCisBoundIndex().toString() + "," +"\n",
                                cadList2.indexOf(cad2),  cad2.getTime(), cad2.getState(), cad2.getCisState(), cad2.getActState(), cad2.isCisM(), cad2.isCisX(), cad2.getBead_position()[0], cad2.getBead_position()[1],
                                cad2.getBead_position()[2], cad2.getCadBoundIndex(), cad2.getForce()));
                    }

                    fw5.write(String.format("time %6.4f \n", cadList1.get(0).getTime()));
                    fw5.write(info.toString() + "\n");
                    fw5.write("--------------------------------------------------- \n");
                    fw5.flush();
                }

                // step 6: record the info in the first second
                if (cadList1.get(0).getTime() < 1 && i % 10000 == 0) {
                    System.out.println("-------------------------------------------------");
                    System.out.println("0.1s record!");
                    double[] count1 = count(cadList1, nameofFile2);
                    System.out.println("the proportion for each state for layer one within one second is: " + Arrays.toString(count1));
                    //System.out.println("the proportion for each state for layer two is: " + Arrays.toString(count2));
                    for (Cadherin cad1 : cadList1) {
                        fw6.write(String.format("layer1-%d," +  "%6.4f," + "%d," + "%d," + "%d," + "%b," + "%b," + "%6.4f," + "%6.4f," + "%6.4f," + "%d," + "%6.4f," + cad1.getCadCisBoundIndex().toString() + "," +"\n",
                                cadList1.indexOf(cad1),  cad1.getTime(), cad1.getState(), cad1.getCisState(), cad1.getActState(), cad1.isCisM(), cad1.isCisX(), cad1.getBead_position()[0], cad1.getBead_position()[1],
                                cad1.getBead_position()[2], cad1.getCadBoundIndex(), cad1.getForce()));
                    }

                    for (Cadherin cad2 : cadList2) {
                        fw6.write(String.format("layer2-%d," + "%6.4f," + "%d," + "%d," + "%d," + "%b," + "%b," + "%6.4f," + "%6.4f," + "%6.4f," + "%d," + "%6.4f," + cad2.getCadCisBoundIndex().toString() + "," + "\n",
                                cadList2.indexOf(cad2), cad2.getTime(), cad2.getState(), cad2.getCisState(), cad2.getActState(), cad2.isCisM(), cad2.isCisX(), cad2.getBead_position()[0], cad2.getBead_position()[1],
                                cad2.getBead_position()[2], cad2.getCadBoundIndex(), cad2.getForce()));
                    }

                    fw7.write(String.format("time %6.4f \n", cadList1.get(0).getTime()));
                    fw7.write(info.toString() + "\n");
                    fw7.write("--------------------------------------------------- \n");
                    fw7.flush();
                }
            }
            fw.close();
            fw5.close();
            fw6.close();
            fw7.close();
        } catch (IOException ioe) {
            System.err.println("IOException: " + ioe.getMessage());
        }
    }

    public static void main(String[] args) {
        // first test can we generate these cadherin objects
        double temperature = 300;// unit: Kelvin
        double domainSize = 1; // unit: micrometer
        double timestep = 10e-6;
        double totalTime = 30; // unit: s
        double zCad_1 = 0;
        double zCad_2 = 0.02; // if we use 20 nm then this value should be 0.02 um
        int nOfCad_1 = 100; // assumption at this stage
        int nOfCad_2 = 100;
        double frictional = 0.1469; // unit: pN*s/um with D = 28e-3 um^2/s

        // generate the Arraylist for Cadherin layer 1 and 2
        ArrayList<Cadherin> Cadherin_1 = new ArrayList<>();
        ArrayList<Cadherin> Cadherin_2 = new ArrayList<>();

        // generate the Arraylist for Actin layer 1 and 2
        ArrayList<Actin> Actin_1 = new ArrayList<>();
        ArrayList<Actin> Actin_2 = new ArrayList<>();

        for (int i = 0; i < nOfCad_1; i++) {
            // initiate the position of each cadherin randomly
            double x = -domainSize / 2 + Math.random() * domainSize;
            double y = -domainSize / 2 + Math.random() * domainSize;
            double[] ini_position = {x, y, zCad_1};
            //double[] ini_position = {0,0,zCad_1}; // start at the origin
            Cadherin_1.add(new Cadherin(ini_position, frictional));
        }

        for (int i = 0; i < nOfCad_2; i++) {
            double x = -domainSize / 2 + Math.random() * domainSize;
            double y = -domainSize / 2 + Math.random() * domainSize;
            double[] ini_position = {x, y, zCad_2};
            Cadherin_2.add(new Cadherin(ini_position, frictional));
        }

        Force_Interaction force_interaction = new Force_Interaction(Cadherin_1, Cadherin_2, Actin_1, Actin_2, timestep);


        // simulate the two layer force interaction
        System.out.println("Start simulation:");
        System.out.printf("Simulate with time step: %8.6f\n", timestep);
        simulate(Cadherin_1, Cadherin_2, totalTime, force_interaction, timestep, "actin100_f50p100_1214#1", temperature, domainSize, frictional);
        System.out.println("Finished!!!");
    }
}
