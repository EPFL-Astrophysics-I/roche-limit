using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RigidSatellite : MonoBehaviour
{
    public bool autoUpdate;

    [Header("Units")]
    [SerializeField] private Units.UnitLength unitLength = Units.UnitLength.EarthRadius;
    [SerializeField] private Units.UnitMass unitMass = Units.UnitMass.EarthMass;
    [SerializeField] private Units.UnitTime unitTime = Units.UnitTime.Day;
    [SerializeField, Min(0)] private float timeScale = 1;
    private double newtonG;

    [Header("Satellite")]
    [SerializeField] private GameObject particlePrefab;
    [SerializeField] private GameObject phantomPrefab;
    [SerializeField] Vector3 initialPosition = 10 * Vector3.left;
    [SerializeField, Min(0)] private int numParticles = 100;
    [SerializeField, Min(0.01f)] private float systemRadius = 10;
    [SerializeField, Min(0.01f)] private float particleRadius = 1;
    private Transform satellite;
    private Transform phantom;
    private double lunarRadius;

    [Header("Primary")]
    [SerializeField] private GameObject primaryPrefab;
    [SerializeField] private float primaryRadius = 1;
    private Transform primary;
    private double primaryMass;

    [Header("Solver")]
    [SerializeField, Min(1)] private int numSubsteps = 1;
    //[SerializeField, Min(0.00001f)] private float epsilon = 0.01f;

    [Header("Controls")]
    [SerializeField] private bool orbiting = true;
    [SerializeField] private bool restoringForces;
    [SerializeField, Min(0)] private double springConstant = 1;
    //[SerializeField, Min(0)] private double springDamping = 0;
    [SerializeField] private bool tidalForces;
    [SerializeField, Min(0)] private double tidalForceFactor = 0;
    [SerializeField] private bool insideRocheLimit;
    [SerializeField] private double kxz = 0;
    [SerializeField] private double ky = 0;

    private double[] x;
    private double[] store;
    private double[] k1;
    //private double[] k2;
    //private double[] k3;
    //private double[] k4;

    // For precise positions of particles in circular orbits
    private double rCM;
    private double yCM;
    private double thetaCM;
    private double omegaCM;
    private double[] thetaRef;
    private double[] rRef;
    private double[] xRef;

    private double G => newtonG;
    private double M => primaryMass;

    public void Clear()
    {
        // Find existing references after exiting Play mode
        satellite = transform.Find("Satellite");
        primary = transform.Find("Primary");
        phantom = transform.Find("Phantom");

        if (satellite)
        {
            // Clear out old particles
            for (int i = satellite.childCount; i > 0; i--)
            {
                DestroyImmediate(satellite.GetChild(0).gameObject);
            }

            DestroyImmediate(satellite.gameObject);
            satellite = null;
        }

        if (primary)
        {
            DestroyImmediate(primary.gameObject);
            primary = null;
        }

        if (phantom)
        {
            // Clear out old phantom particles
            for (int i = phantom.childCount; i > 0; i--)
            {
                DestroyImmediate(phantom.GetChild(0).gameObject);
            }

            DestroyImmediate(phantom.gameObject);
            phantom = null;
        }
    }

    public void Generate()
    {
        Clear();

        if (!satellite)
        {
            satellite = new GameObject("Satellite").transform;
            satellite.parent = transform;
        }

        if (!particlePrefab)
        {
            Debug.LogWarning("No Piece Prefab assigned in RigidSatellite");
            return;
        }

        // Create new particles
        for (int i = 0; i < numParticles; i++)
        {
            Transform particle = Instantiate(particlePrefab, satellite).transform;
            particle.localScale = 2 * particleRadius * Vector3.one;
        }

        // Distribute evenly around a sphere
        Vector3 positionCM = Vector3.zero;

        float turnFraction = 0.5f * (1 + Mathf.Sqrt(5));  // golden ratio

        // Evenly space N bodies around a sphere
        for (int i = 0; i < numParticles; i++)
        {
            Vector3 position = Vector3.zero;
            if (numParticles > 1)
            {
                float t = i / (numParticles - 1f);
                float inclination = Mathf.Acos(1 - 2 * t);
                float azimuth = 2 * Mathf.PI * turnFraction * i;

                float positionX = Mathf.Sin(inclination) * Mathf.Cos(azimuth);
                float positionY = Mathf.Sin(inclination) * Mathf.Sin(azimuth);
                float positionZ = Mathf.Cos(inclination);
                position = systemRadius * new Vector3(positionX, positionY, positionZ);
            }

            satellite.GetChild(i).position = position;
            positionCM += position;
        }

        positionCM /= numParticles;

        // Work in the CM frame (i.e. shift the system to the origin)
        for (int i = 0; i < numParticles; i++)
        {
            satellite.GetChild(i).position += initialPosition - positionCM;
        }

        if (primaryPrefab)
        {
            primary = Instantiate(primaryPrefab, transform).transform;
            primary.localScale = 2 * primaryRadius * Vector3.one;
            primary.name = "Primary";
        }

        if (phantomPrefab)
        {
            phantom = Instantiate(phantomPrefab, transform).transform;
            phantom.position = initialPosition;
            phantom.localScale = Vector3.one;
            phantom.name = "Phantom";

            for (int i = 0; i < numParticles; i++)
            {
                Transform phantomParticle = Instantiate(phantomPrefab, phantom).transform;
                phantomParticle.localScale = particleRadius * Vector3.one;
                phantomParticle.position = satellite.GetChild(i).position;
            }
        }
    }

    private void Start()
    {
        // Find existing references when entering Play mode
        satellite = transform.Find("Satellite");
        primary = transform.Find("Primary");
        phantom = transform.Find("Phantom");

        // Compute gravitational constant
        newtonG = Units.NewtonG(unitTime, unitLength, unitMass);
        //Debug.Log("G = " + G);
        primaryMass = Units.EarthMass(unitMass);
        //Debug.Log("M = " + M);
        lunarRadius = Units.LunarRadius(unitLength);
        //Debug.Log("R = " + lunarRadius);

        if (!satellite)
        {
            return;
        }

        int numBodies = satellite.childCount;

        // For true equations of motion of particles
        int numEquations = 6 * numBodies;
        x = new double[numEquations];
        store = new double[numEquations];
        k1 = new double[numEquations];
        //k2 = new double[numEquations];
        //k3 = new double[numEquations];
        //k4 = new double[numEquations];

        // For reference phantom particles in uniform circular motion
        thetaRef = new double[numBodies];
        rRef = new double[numBodies];
        xRef = new double[3 * numBodies];

        // Compute satellite CM angular speed
        thetaCM = Mathf.Atan2(initialPosition.z, initialPosition.x);
        rCM = Mathf.Sqrt(initialPosition.x * initialPosition.x + initialPosition.z + initialPosition.z);
        yCM = initialPosition.y;
        omegaCM = System.Math.Sqrt(G * M / rCM / rCM / rCM);
        Debug.Log("P = " + (2 * Mathf.PI / omegaCM));

        for (int i = 0; i < numBodies; i++)
        {
            Transform particle = satellite.GetChild(i);

            // Initial position
            x[i * 6 + 0] = particle.position.x;
            x[i * 6 + 1] = particle.position.y;
            x[i * 6 + 2] = particle.position.z;

            // Initial velocity
            double r = System.Math.Sqrt(x[i * 6 + 0] * x[i * 6 + 0] + x[i * 6 + 2] * x[i * 6 + 2]);
            double speed = omegaCM * r;
            x[i * 6 + 3] = -speed * x[i * 6 + 2] / r;
            x[i * 6 + 4] = 0;
            x[i * 6 + 5] = speed * x[i * 6 + 0] / r;

            thetaRef[i] = System.Math.Atan2(x[i * 6 + 2], x[i * 6 + 0]);
            rRef[i] = r;

            // Save initial position and velocity
            xRef[i * 3 + 0] = x[i * 6 + 0];
            xRef[i * 3 + 1] = x[i * 6 + 1];
            xRef[i * 3 + 2] = x[i * 6 + 2];
        }

        // Compute accelerations for a valid first Verlet step
        RateOfChange(x, store);
    }

    private void FixedUpdate()
    {
        if (!satellite)
        {
            return;
        }

        // Move the reference
        if (orbiting)
        {
            double deltaTheta = omegaCM * timeScale * Time.fixedDeltaTime;
            thetaCM += deltaTheta;
            if (thetaCM >= 2 * System.Math.PI)
            {
                thetaCM -= 2 * System.Math.PI;
            }
            //Debug.Log(thetaCM * 180 / Mathf.PI);
            for (int i = 0; i < satellite.childCount; i++)
            {
                thetaRef[i] += deltaTheta;
                xRef[i * 3 + 0] = rRef[i] * System.Math.Cos(thetaRef[i]);
                xRef[i * 3 + 2] = rRef[i] * System.Math.Sin(thetaRef[i]);
            }
        }

        float substep = timeScale * Time.fixedDeltaTime / numSubsteps;
        for (int i = 0; i < numSubsteps; i++)
        {
            VerletStep(substep);
        }

        // Update game object positions
        for (int i = 0; i < satellite.childCount; i++)
        {
            satellite.GetChild(i).position = new Vector3((float)x[i * 6 + 0], (float)x[i * 6 + 1], (float)x[i * 6 + 2]);
        }

        if (phantom)
        {
            double xCM = rCM * System.Math.Cos(thetaCM);
            double zCM = rCM * System.Math.Sin(thetaCM);
            phantom.position = new Vector3((float)xCM, (float)yCM, (float)zCM);

            for (int i = 0; i < phantom.childCount; i++)
            {
                phantom.GetChild(i).position = new Vector3((float)xRef[i * 3 + 0], (float)xRef[i * 3 + 1], (float)xRef[i * 3 + 2]);
            }
        }
    }

    private void VerletStep(float deltaTime)
    {
        int numBodies = satellite.childCount;
        for (int i = 0; i < numBodies; i++)
        {
            // Update positions
            x[i * 6 + 0] += deltaTime * (x[i * 6 + 3] + 0.5f * deltaTime * store[i * 6 + 3]);
            x[i * 6 + 1] += deltaTime * (x[i * 6 + 4] + 0.5f * deltaTime * store[i * 6 + 4]);
            x[i * 6 + 2] += deltaTime * (x[i * 6 + 5] + 0.5f * deltaTime * store[i * 6 + 5]);
        }

        // Compute accelerations based on updated positions
        RateOfChange(x, k1);

        for (int i = 0; i < numBodies; i++)
        {
            // Update velocities
            x[i * 6 + 3] += 0.5f * deltaTime * (store[i * 6 + 3] + k1[i * 6 + 3]);
            x[i * 6 + 4] += 0.5f * deltaTime * (store[i * 6 + 4] + k1[i * 6 + 4]);
            x[i * 6 + 5] += 0.5f * deltaTime * (store[i * 6 + 5] + k1[i * 6 + 5]);

            // Keep previous acceleration values
            store[i * 6 + 3] = k1[i * 6 + 3];
            store[i * 6 + 4] = k1[i * 6 + 4];
            store[i * 6 + 5] = k1[i * 6 + 5];
        }
    }

    private void RateOfChange(double[] x, double[] xdot)
    {
        int numBodies = satellite.childCount;

        // TODO compute vector to satellite CM position outside the loop

        // Set xdot and vdot
        for (int i = 0; i < numBodies; i++)
        {
            xdot[i * 6 + 0] = x[i * 6 + 3];
            xdot[i * 6 + 1] = x[i * 6 + 4];
            xdot[i * 6 + 2] = x[i * 6 + 5];
            // Start with zero acceleration; add forces one by one as required
            xdot[i * 6 + 3] = 0;
            xdot[i * 6 + 4] = 0;
            xdot[i * 6 + 5] = 0;

            if (insideRocheLimit)
            {
                // Particles only feel the gravitational force of the primary (in the x-z plane),
                // plus a restoring force in the y-direction
                double rx = x[i * 6 + 0];
                double rz = x[i * 6 + 2];
                double r2 = rx * rx + rz * rz;
                double r = System.Math.Sqrt(r2);
                // Acceleration due to uniform circular motion
                if (r > 0)
                {
                    xdot[i * 6 + 3] += -G * M / r2 * rx / r - kxz * (1 - rCM / r) * rx;
                    xdot[i * 6 + 4] += -ky * x[i * 6 + 1] - 2 * System.Math.Sqrt(ky) * xdot[i * 6 + 1];
                    xdot[i * 6 + 5] += -G * M / r2 * rz / r - kxz * (1 - rCM / r) * rz;
                }
            }
            else
            {
                // Satellite orbits in the x-z plane about the origin
                if (orbiting)
                {
                    // Projected velocity in the x-z plane
                    double vx = xdot[i * 6 + 0];
                    double vz = xdot[i * 6 + 2];
                    double v2 = vx * vx + vz * vz;
                    // Projected position in the x-z plane
                    double rx = x[i * 6 + 0];
                    double rz = x[i * 6 + 2];
                    double r2 = rx * rx + rz * rz;
                    // Acceleration due to uniform circular motion
                    if (r2 > 0)
                    {
                        xdot[i * 6 + 3] += -v2 * rx / r2;
                        xdot[i * 6 + 5] += -v2 * rz / r2;
                    }
                }

                // Restoring forces bringing particles back to (true) reference positions
                if (restoringForces)
                {
                    double dx = x[i * 6 + 0] - xRef[i * 3 + 0];
                    double dy = x[i * 6 + 1] - xRef[i * 3 + 1];
                    double dz = x[i * 6 + 2] - xRef[i * 3 + 2];

                    // True orbital reference speed
                    //double speed = omegaCM * rRef[i];
                    double vx = -omegaCM * x[i * 6 + 2];
                    double vy = 0;
                    double vz = omegaCM * x[i * 6 + 0];

                    double damping = 2 * System.Math.Sqrt(springConstant);
                    xdot[i * 6 + 3] += -springConstant * dx - damping * (xdot[i * 6 + 0] - vx);
                    xdot[i * 6 + 4] += -springConstant * dy - damping * (xdot[i * 6 + 1] - vy);
                    xdot[i * 6 + 5] += -springConstant * dz - damping * (xdot[i * 6 + 2] - vz);
                }

                if (tidalForces)
                {

                    double factor = tidalForceFactor * G * M * lunarRadius;

                    // Satellite CM position
                    // TODO what to do with y position ??
                    double rCMx = rCM * System.Math.Cos(thetaCM);
                    double rCMy = 0;
                    double rCMz = rCM * System.Math.Sin(thetaCM);
                    double rCM2 = rCMx * rCMx + rCMz * rCMz;
                    double rCMmag = System.Math.Sqrt(rCM2);

                    // Particle position in the satellite CM frame
                    double bx = x[i * 6 + 0] - rCMx;
                    double by = x[i * 6 + 1] - rCMy;
                    double bz = x[i * 6 + 2] - rCMz;
                    double b2 = bx * bx + by * by + bz * bz;
                    double b = System.Math.Sqrt(b2);

                    double cosTheta = (bx * rCMx + by * rCMy + bz * rCMz) / b / rCMmag;
                    double sinTheta = System.Math.Sqrt(1 - cosTheta * cosTheta);

                    // First term (i-hat = -r-hat)
                    xdot[i * 6 + 3] += factor * 2 * cosTheta * rCMx / rCM2 / rCMmag;
                    xdot[i * 6 + 4] += factor * 2 * cosTheta * rCMy / rCM2 / rCMmag;
                    xdot[i * 6 + 5] += factor * 2 * cosTheta * rCMx / rCM2 / rCMmag;
                }
            }

            //if (!insideRocheLimit)
            //{
            //    double dx = x[i * 6 + 0] - x0[i * 6 + 0];
            //    double dy = x[i * 6 + 1] - x0[i * 6 + 1];
            //    double dz = x[i * 6 + 2] - x0[i * 6 + 2];

            //    xdot[i * 6 + 3] += -springConstant * dx - springDamping * xdot[i * 6 + 0];
            //    xdot[i * 6 + 4] += -springConstant * dy - springDamping * xdot[i * 6 + 1];
            //    xdot[i * 6 + 5] += -springConstant * dz - springDamping * xdot[i * 6 + 2];

            //    if (tidalForces && primary != null)
            //    {
            //        double factor = G * M * systemRadius;

            //        // Vector to primary
            //        double rx = primary.position.x - transform.position.x;
            //        double ry = primary.position.y - transform.position.y;
            //        double rz = primary.position.z - transform.position.z;
            //        double r2 = rx * rx + ry * ry + rz * rz;
            //        double r = System.Math.Sqrt(r2);

            //        // Vector to body
            //        double bx = x[i * 6 + 0];
            //        double by = x[i * 6 + 1];
            //        double bz = x[i * 6 + 2];
            //        double b2 = bx * bx + by * by + bz * bz;
            //        double b = System.Math.Sqrt(b2);

            //        double cosTheta = (bx * rx + by * ry + bz * rz) / b / r;
            //        double sinTheta = System.Math.Sqrt(1 - cosTheta * cosTheta);

            //        // First term (r-hat = i-hat)
            //        xdot[i * 6 + 3] += factor * 2 * cosTheta * rx / r2 / r;
            //        xdot[i * 6 + 4] += factor * 2 * cosTheta * ry / r2 / r;
            //        xdot[i * 6 + 5] += factor * 2 * cosTheta * rz / r2 / r;

            //        // Second term
            //        if (cosTheta != 0)
            //        {
            //            // Find local j-hat vector
            //            double tx = (ry * bz - rz * by) / b / r;
            //            double ty = (rz * bx - rx * bz) / b / r;
            //            double tz = (rx * by - ry * bx) / b / r;

            //            double jx = (ty * rz - tz * ry) / r;
            //            double jy = (tz * rx - tx * rz) / r;
            //            double jz = (tx * ry - ty * rx) / r;

            //            xdot[i * 6 + 3] -= factor * sinTheta * jx / r2 / r;
            //            xdot[i * 6 + 4] -= factor * sinTheta * jy / r2 / r;
            //            xdot[i * 6 + 5] -= factor * sinTheta * jz / r2 / r;
            //        }
            //    }
            //}
            //else
            //{
            //    double dx = x[i * 6 + 0] - (primary.position.x - transform.position.x);
            //    double dy = x[i * 6 + 1] - (primary.position.y - transform.position.y);
            //    double dz = x[i * 6 + 2] - (primary.position.z - transform.position.z);
            //    double dr2 = dx * dx + dy * dy + dz * dz + epsilon * epsilon;
            //    double dr = System.Math.Sqrt(dr2);
            //    xdot[i * 6 + 3] += -G * M / dr2 * dx / dr;
            //    xdot[i * 6 + 4] += -G * M / dr2 * dy / dr;
            //    xdot[i * 6 + 5] += -G * M / dr2 * dz / dr;
            //}
        }
    }
}
