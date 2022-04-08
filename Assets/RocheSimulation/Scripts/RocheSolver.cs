using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RocheSolver
{
    public int numBodies;
    private double newtonG;
    public double[] x;  // positions and velocities
    private double mass;  // universal mass
    private double epsilon;  // softening

    // Integration arrays
    private int numEquations;
    private double[] store;
    private double[] k1;
    private double[] k2;
    private double[] k3;
    private double[] k4;

    // Constructor
    public RocheSolver(int numBodies, double newtonG, double mass, double epsilon)
    {
        this.numBodies = numBodies;
        this.newtonG = newtonG;
        this.mass = mass;
        this.epsilon = epsilon;

        numEquations = 6 * numBodies;
        x = new double[numEquations];
        store = new double[numEquations];
        k1 = new double[numEquations];
        k2 = new double[numEquations];
        k3 = new double[numEquations];
        k4 = new double[numEquations];
    }

    public void Initialize(float radialMean, float radialSigma, Vector3 systemPosition)
    {
        Vector3 positionCM = Vector3.zero;

        float turnFraction = 0.5f * (1 + Mathf.Sqrt(5));  // golden ratio

        // Evenly space N bodies around a sphere
        for (int i = 0; i < numBodies; i++)
        {
            float t = i / (numBodies - 1f);
            float inclination = Mathf.Acos(1 - 2 * t);
            float azimuth = 2 * Mathf.PI * turnFraction * i;

            float positionX = Mathf.Sin(inclination) * Mathf.Cos(azimuth);
            float positionY = Mathf.Sin(inclination) * Mathf.Sin(azimuth);
            float positionZ = Mathf.Cos(inclination);

            float radius = Mathf.Abs(Utils.Random.NormalValue(radialMean, radialSigma));
            //Vector3 position = radius * Random.onUnitSphere;
            //Vector3 position = radius * new Vector3(positionX, positionY, positionZ);
            //positions.Add(position);
            x[i * 6 + 0] = radius * positionX;
            x[i * 6 + 1] = radius * positionY;
            x[i * 6 + 2] = radius * positionZ;
            positionCM += radius * new Vector3(positionX, positionY, positionZ);
        }

        positionCM /= numBodies;

        // Work in the CM frame (i.e. shift the system to the origin)
        for (int i = 0; i < numBodies; i++)
        {
            //positions[i] -= positionCM;
            x[i * 6 + 0] += systemPosition.x - positionCM.x;
            x[i * 6 + 1] += systemPosition.y - positionCM.y;
            x[i * 6 + 2] += systemPosition.z - positionCM.z;
        }
    }

    public void StepForward(double deltaTime, NBodySolver.IntegrationMethod method)
    {

    }
}
