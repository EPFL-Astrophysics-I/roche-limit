using UnityEngine;

[RequireComponent(typeof(MoonRingPrefabs))]
public class MoonRing : Simulation
{
    private MoonRingPrefabs prefabs;

    [Header("Units")]
    [SerializeField] private Units.UnitLength unitLength = Units.UnitLength.EarthRadius;
    [SerializeField] private Units.UnitMass unitMass = Units.UnitMass.EarthMass;
    [SerializeField] private Units.UnitTime unitTime = Units.UnitTime.Day;
    [SerializeField, Min(0)] private float timeScale = 1;
    private double newtonG;

    [Header("Earth")]
    [SerializeField] private float earthRadius = 1;
    [SerializeField] private bool earthIsRotating = false;
    private double earthMass;
    private float earthRotationFrequency;

    [Header("Moon")]
    public float lunarDistance = 10;
    [SerializeField, Min(0)] private int numParticles = 100;
    [SerializeField, Min(0.01f)] private float particleRadius = 1;
    private double lunarRadius;
    private Vector3 initialMoonPosition;

    [Header("Solver")]
    [SerializeField, Min(1)] private int numSubsteps = 1;
    //[SerializeField, Min(0.00001f)] private float epsilon = 0.01f;

    [Header("Controls")]
    public bool moonIsSolid = true;
    //[SerializeField] private bool restoringForces;
    //[SerializeField, Min(0)] private double springConstant = 1;
    [SerializeField] private double kxz = 0;
    [SerializeField] private double ky = 0;

    private double[] x;
    private double[] store;
    private double[] k1;
    //private double[] k2;
    //private double[] k3;
    //private double[] k4;

    // For precise positions of particles in circular orbits
    //private double rCM;
    //private double yCM;
    private double thetaCM;
    private double omegaCM;
    //private double[] thetaRef;
    //private double[] rRef;
    //private double[] xRef;

    private double G => newtonG;
    private double M => earthMass;
    private double R => lunarRadius;

    private void Awake()
    {
        // Set the gravitational constant for these units
        newtonG = Units.NewtonG(unitTime, unitLength, unitMass);        
        earthMass = Units.EarthMass(unitMass);
        lunarRadius = Units.LunarRadius(unitLength);
        //Debug.Log("G = " + G);
        //Debug.Log("M = " + M);
        //Debug.Log("R = " + R);

        // Create all objects with assigned prefabs
        if (transform.TryGetComponent(out prefabs))
        {
            prefabs.InstantiatePrefabs(earthRadius, (float)(earthRadius * R), lunarDistance);
        }

        // Save initial moon position for resetting
        if (prefabs.solidMoon)
        {
            initialMoonPosition = prefabs.solidMoon.position;
        }

        if (prefabs.earth)
        {
            prefabs.earth.localRotation = Quaternion.Euler(115 * Vector3.up);
        }
    }

    private void ResetSolidMoon()
    {
        // Compute the moon's (constant) angular speed
        if (prefabs.solidMoon)
        {
            Vector3 position = initialMoonPosition;
            thetaCM = Mathf.Atan2(position.z, position.x);
            float r = Mathf.Sqrt(position.x * position.x + position.z + position.z);
            //yCM = position.y;
            omegaCM = System.Math.Sqrt(G * M / r / r / r);

            float orbitalPeriod = 2 * Mathf.PI / (float)omegaCM;
            Debug.Log("P = " + orbitalPeriod);
            earthRotationFrequency = 360 * 27 / orbitalPeriod;

            prefabs.solidMoon.gameObject.SetActive(true);
        }
    }

    private void ClearMoonParticles()
    {
        if (prefabs.moonParticles)
        {
            Destroy(prefabs.moonParticles.gameObject);
            prefabs.moonParticles = null;
        }
    }

    private void Start()
    {
        ResetSolidMoon();
        ClearMoonParticles();
        moonIsSolid = true;
    }

    private void FixedUpdate()
    {
        if (prefabs.earth && earthIsRotating)
        {
            //float deltaTheta = 0.1f;
            float deltaTheta = earthRotationFrequency * timeScale * Time.fixedDeltaTime;
            prefabs.earth.RotateAround(Vector3.zero, Vector3.down, deltaTheta);
        }

        if (prefabs.solidMoon && moonIsSolid)
        {
            // Update the angular position of the solid moon
            double deltaTheta = omegaCM * timeScale * Time.fixedDeltaTime;
            thetaCM += deltaTheta;
            if (thetaCM >= 2 * System.Math.PI)
            {
                thetaCM -= 2 * System.Math.PI;
            }

            Vector3 direction = new Vector3((float)System.Math.Cos(thetaCM), 0, (float)System.Math.Sin(thetaCM));
            prefabs.solidMoon.position = lunarDistance * direction;
            prefabs.solidMoon.localRotation = Quaternion.Euler((float)thetaCM * 180f / Mathf.PI * Vector3.down);
        }

        if (prefabs.moonParticles)
        {
            float substep = timeScale * Time.fixedDeltaTime / numSubsteps;
            for (int i = 0; i < numSubsteps; i++)
            {
                VerletStep(substep);
            }

            // Update game object positions
            for (int i = 0; i < prefabs.moonParticles.childCount; i++)
            {
                prefabs.moonParticles.GetChild(i).position = new Vector3((float)x[i * 6 + 0], (float)x[i * 6 + 1], (float)x[i * 6 + 2]);
            }
        }

        //if (moonIsSolid)
        //{
        //    // Update the angular position of the solid moon
        //    double deltaTheta = omegaCM * timeScale * Time.fixedDeltaTime;
        //    thetaCM += deltaTheta;
        //    if (thetaCM >= 2 * System.Math.PI)
        //    {
        //        thetaCM -= 2 * System.Math.PI;
        //    }

        //    // Compute exact 3D particle positions
        //    for (int i = 0; i < thetaRef.Length; i++)
        //    {
        //        thetaRef[i] += deltaTheta;
        //        xRef[i * 3 + 0] = rRef[i] * System.Math.Cos(thetaRef[i]);
        //        xRef[i * 3 + 2] = rRef[i] * System.Math.Sin(thetaRef[i]);
        //    }

        //    // Update the phantom particle transforms
        //    if (prefabs.phantom)
        //    {
        //        double xCM = rCM * System.Math.Cos(thetaCM);
        //        double zCM = rCM * System.Math.Sin(thetaCM);
        //        prefabs.phantom.position = new Vector3((float)xCM, (float)yCM, (float)zCM);

        //        for (int i = 0; i < prefabs.phantom.childCount; i++)
        //        {
        //            prefabs.phantom.GetChild(i).position = new Vector3((float)xRef[i * 3 + 0], (float)xRef[i * 3 + 1], (float)xRef[i * 3 + 2]);
        //        }
        //    }
        //}
    }

    private void VerletStep(float deltaTime)
    {
        for (int i = 0; i < prefabs.numMoonParticles; i++)
        {
            // Update positions
            x[i * 6 + 0] += deltaTime * (x[i * 6 + 3] + 0.5f * deltaTime * store[i * 6 + 3]);
            x[i * 6 + 1] += deltaTime * (x[i * 6 + 4] + 0.5f * deltaTime * store[i * 6 + 4]);
            x[i * 6 + 2] += deltaTime * (x[i * 6 + 5] + 0.5f * deltaTime * store[i * 6 + 5]);
        }

        // Compute accelerations based on updated positions
        RateOfChange(x, k1);

        for (int i = 0; i < prefabs.numMoonParticles; i++)
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
        // TODO compute vector to satellite CM position outside the loop

        // Set xdot and vdot
        for (int i = 0; i < prefabs.numMoonParticles; i++)
        {
            xdot[i * 6 + 0] = x[i * 6 + 3];
            xdot[i * 6 + 1] = x[i * 6 + 4];
            xdot[i * 6 + 2] = x[i * 6 + 5];
            // Start with zero acceleration; add forces one by one as required
            xdot[i * 6 + 3] = 0;
            xdot[i * 6 + 4] = 0;
            xdot[i * 6 + 5] = 0;

            //if (moonIsSolid)
            //{
            //    // Moon orbits in the x-z plane about the origin
            //    double vx = xdot[i * 6 + 0];
            //    double vz = xdot[i * 6 + 2];
            //    double v2 = vx * vx + vz * vz;
            //    // Projected position in the x-z plane
            //    double rx = x[i * 6 + 0];
            //    double rz = x[i * 6 + 2];
            //    double r2 = rx * rx + rz * rz;
            //    // Acceleration due to uniform circular motion
            //    if (r2 > 0)
            //    {
            //        xdot[i * 6 + 3] += -v2 * rx / r2;
            //        xdot[i * 6 + 5] += -v2 * rz / r2;
            //    }

            //    // Restoring forces bringing particles back to (true) reference positions
            //    if (restoringForces)
            //    {
            //        double dx = x[i * 6 + 0] - xRef[i * 3 + 0];
            //        double dy = x[i * 6 + 1] - xRef[i * 3 + 1];
            //        double dz = x[i * 6 + 2] - xRef[i * 3 + 2];

            //        // True orbital reference speed
            //        //double speed = omegaCM * rRef[i];
            //        double ux = -omegaCM * x[i * 6 + 2];
            //        double uy = 0;
            //        double uz = omegaCM * x[i * 6 + 0];

            //        double damping = 2 * System.Math.Sqrt(springConstant);
            //        xdot[i * 6 + 3] += -springConstant * dx - damping * (xdot[i * 6 + 0] - ux);
            //        xdot[i * 6 + 4] += -springConstant * dy - damping * (xdot[i * 6 + 1] - uy);
            //        xdot[i * 6 + 5] += -springConstant * dz - damping * (xdot[i * 6 + 2] - uz);
            //    }
            //}
            //else
            //{
            //    // Particles only feel the gravitational force of the primary (in the x-z plane),
            //    // plus a restoring force in the y-direction
            //    double rx = x[i * 6 + 0];
            //    double rz = x[i * 6 + 2];
            //    double r2 = rx * rx + rz * rz;
            //    double r = System.Math.Sqrt(r2);
            //    // Acceleration due to uniform circular motion
            //    if (r > 0)
            //    {
            //        xdot[i * 6 + 3] += -G * M / r2 * rx / r - kxz * (1 - rCM / r) * rx;
            //        xdot[i * 6 + 4] += -ky * x[i * 6 + 1] - 2 * System.Math.Sqrt(ky) * xdot[i * 6 + 1];
            //        xdot[i * 6 + 5] += -G * M / r2 * rz / r - kxz * (1 - rCM / r) * rz;
            //    }
            //}

            // Particles only feel the gravitational force of the primary (in the x-z plane),
            // plus a restoring force in the y-direction
            double rx = x[i * 6 + 0];
            double rz = x[i * 6 + 2];
            double r2 = rx * rx + rz * rz;
            double r = System.Math.Sqrt(r2);
            if (r > 0)
            {
                // Acceleration due to Earth's gravity
                xdot[i * 6 + 3] += -G * M / r2 * rx / r - kxz * (1 - lunarDistance / r) * rx;
                xdot[i * 6 + 4] += -ky * x[i * 6 + 1] - 2 * System.Math.Sqrt(ky) * xdot[i * 6 + 1];
                xdot[i * 6 + 5] += -G * M / r2 * rz / r - kxz * (1 - lunarDistance / r) * rz;
            }
        }
    }

    public void BreakUpMoon()
    {
        if (!prefabs.solidMoon)
        {
            return;
        }

        prefabs.solidMoon.gameObject.SetActive(false);

        if (prefabs.rocheLimitLR)
        {
            prefabs.rocheLimitLR.gameObject.SetActive(false);
        }

        prefabs.GenerateMoonParticles(particleRadius, (float)(earthRadius * R), prefabs.solidMoon.position);

        // Reference phantom particles in uniform circular motion
        //thetaRef = new double[numParticles];
        //rRef = new double[numParticles];
        //xRef = new double[3 * numParticles];

        //if (!prefabs.solidMoon)
        //{
        //    return;
        //}

        // For true equations of motion of particles
        int numEquations = 6 * prefabs.numMoonParticles;
        x = new double[numEquations];
        store = new double[numEquations];
        k1 = new double[numEquations];
        //k2 = new double[numEquations];
        //k3 = new double[numEquations];
        //k4 = new double[numEquations];

        // Compute satellite CM angular speed
        //Vector3 initialPosition = lunarDistance * Vector3.left;
        //thetaCM = Mathf.Atan2(initialPosition.z, initialPosition.x);
        //rCM = Mathf.Sqrt(initialPosition.x * initialPosition.x + initialPosition.z + initialPosition.z);
        //yCM = initialPosition.y;
        //omegaCM = System.Math.Sqrt(G * M / rCM / rCM / rCM);
        //Debug.Log("P = " + (2 * Mathf.PI / omegaCM));

        for (int i = 0; i < prefabs.numMoonParticles; i++)
        {
            Transform particle = prefabs.moonParticles.GetChild(i);

            // Initial position
            x[i * 6 + 0] = particle.position.x;
            x[i * 6 + 1] = particle.position.y;
            x[i * 6 + 2] = particle.position.z;

            // Initial velocity
            double r = System.Math.Sqrt(x[i * 6 + 0] * x[i * 6 + 0] + x[i * 6 + 2] * x[i * 6 + 2]);
            double speed = omegaCM * r * (2 - r / lunarDistance);
            x[i * 6 + 3] = -speed * x[i * 6 + 2] / r;
            x[i * 6 + 4] = 0;
            x[i * 6 + 5] = speed * x[i * 6 + 0] / r;

            //thetaRef[i] = System.Math.Atan2(x[i * 6 + 2], x[i * 6 + 0]);
            //rRef[i] = r;

            // Save initial position and velocity
            //xRef[i * 3 + 0] = x[i * 6 + 0];
            //xRef[i * 3 + 1] = x[i * 6 + 1];
            //xRef[i * 3 + 2] = x[i * 6 + 2];
        }

        // Compute accelerations for a valid first Verlet step
        RateOfChange(x, store);
    }

    public override void Reset()
    {
        ResetSolidMoon();
        ClearMoonParticles();

        if (prefabs.rocheLimitLR)
        {
            prefabs.rocheLimitLR.gameObject.SetActive(true);
        }

        moonIsSolid = true;
    }

    public void SetEarthRotating(bool value)
    {
        earthIsRotating = value;
    }
}
