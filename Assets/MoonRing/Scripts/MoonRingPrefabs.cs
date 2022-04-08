using UnityEngine;

public class MoonRingPrefabs : MonoBehaviour
{
    [Header("Earth")]
    [SerializeField] private GameObject earthPrefab;
    [HideInInspector] public Transform earth;
    [HideInInspector] public float earthRadius;

    [Header("Moon")]
    [SerializeField] private GameObject solidMoonPrefab;
    [SerializeField] private GameObject moonParticlePrefab;
    //[SerializeField] private GameObject phantomPrefab;
    [HideInInspector] public Transform solidMoon;
    [HideInInspector] public Transform moonParticles;
    //[HideInInspector] public Transform phantom;
    [HideInInspector] public int numMoonParticles;

    [Header("Roche Limit")]
    [SerializeField] private GameObject rocheLimitPrefab;
    [HideInInspector] public LineRenderer rocheLimitLR;

    [Header("Lights")]
    [SerializeField] private GameObject[] lightPrefabs;
    [HideInInspector] public Transform[] lights;

    // Saved values for resetting
    //private int numParticles;
    //private float moonRadius;
    //private Vector3 moonPosition;

    public void GenerateEarth(float radius)
    {
        if (earthPrefab)
        {
            earth = Instantiate(earthPrefab, transform).transform;
            earth.position = Vector3.zero;
            earth.localScale = 2 * radius * Vector3.one;
            earth.name = "Earth";
            earthRadius = radius;
        }
        else
        {
            Debug.LogWarning("Cannot generate Earth: no prefab assigned.");
        }
    }

    public void GenerateSolidMoon(float radius, float distance)
    {
        if (solidMoonPrefab)
        {
            solidMoon = Instantiate(solidMoonPrefab, transform).transform;
            solidMoon.position = distance * Vector3.left;
            solidMoon.localScale = 2 * radius * Vector3.one;
            solidMoon.name = "Solid Moon";
        }
        else
        {
            Debug.LogWarning("Cannot generate Solid Moon: no prefab assigned.");
        }
    }

    public void GenerateMoonParticles(float particleRadius, float moonRadius, Vector3 moonPosition)
    {
        if (!moonParticlePrefab)
        {
            return;
        }

        // Compute the number of layers needed
        int numLayers = Mathf.FloorToInt(0.5f * (moonRadius / particleRadius + 1));
        Debug.Log("N layers = " + numLayers);

        // Compute the total number of particles
        int numParticles = 1;
        for (int i = 2; i < numLayers + 1; i++)
        {
            numParticles += Mathf.FloorToInt(0.74f * (Mathf.Pow(2 * i - 1, 3) - Mathf.Pow(2 * i - 3, 3)));
        }
        Debug.Log("N part " + numParticles);

        // Create a container for the particles
        moonParticles = new GameObject("Moon Particles").transform;
        moonParticles.SetParent(transform);
        moonParticles.position = Vector3.zero;

        // Create the particles
        for (int i = 0; i < numParticles; i++)
        {
            Transform particle = Instantiate(moonParticlePrefab, moonParticles).transform;
            particle.localPosition = Vector3.zero;
            particle.localScale = 2 * particleRadius * Vector3.one;
        }

        // Distribute the particles roughly evenly throughout a sphere
        Vector3 positionCM = Vector3.zero;

        float turnFraction = 0.5f * (1 + Mathf.Sqrt(5));  // golden ratio

        int particleIndex = 1;
        for (int j = 2; j < numLayers + 1; j++)
        {
            int numParticlesInLayer = Mathf.FloorToInt(0.74f * (Mathf.Pow(2 * j - 1, 3) - Mathf.Pow(2 * j - 3, 3)));
            float layerRadius = 2 * (j - 1) * particleRadius;

            // Evenly space N bodies around a sphere
            for (int i = 0; i < numParticlesInLayer; i++)
            {
                Vector3 particlePosition = Vector3.zero;
                if (numParticlesInLayer > 1)
                {
                    float t = i / (numParticlesInLayer - 1f);
                    float inclination = Mathf.Acos(1 - 2 * t);
                    float azimuth = 2 * Mathf.PI * turnFraction * i;

                    float positionX = Mathf.Sin(inclination) * Mathf.Cos(azimuth);
                    float positionY = Mathf.Sin(inclination) * Mathf.Sin(azimuth);
                    float positionZ = Mathf.Cos(inclination);
                    particlePosition = Random.Range(0.9f, 1.1f) * layerRadius * new Vector3(positionX, positionY, positionZ);
                }

                moonParticles.GetChild(particleIndex).position = particlePosition;
                positionCM += particlePosition;
                particleIndex++;
            }
        }

        positionCM /= numParticles;

        // Work in the CM frame (i.e. shift the system to the origin)
        for (int i = 0; i < numParticles; i++)
        {
            moonParticles.GetChild(i).position += moonPosition - positionCM;
        }

        //if (phantomPrefab)
        //{
        //    phantom = Instantiate(phantomPrefab, transform).transform;
        //    phantom.position = moonPosition;
        //    phantom.localScale = Vector3.one;
        //    phantom.name = "Phantom";

        //    for (int i = 0; i < numParticles; i++)
        //    {
        //        Transform phantomParticle = Instantiate(phantomPrefab, phantom).transform;
        //        phantomParticle.localScale = particleRadius * Vector3.one;
        //        phantomParticle.position = moonParticles.GetChild(i).position;
        //    }
        //}

        numMoonParticles = numParticles;

        //this.numParticles = numParticles;
        //this.moonRadius = moonRadius;
        //this.moonPosition = moonPosition;
    }

    public void DrawRocheLimit(float distance, int numSamples = 1000)
    {
        if (!rocheLimitLR)
        {
            if (rocheLimitPrefab)
            {
                rocheLimitLR = Instantiate(rocheLimitPrefab, transform).GetComponent<LineRenderer>();
                rocheLimitLR.transform.position = Vector3.zero;
                rocheLimitLR.name = "Roche Limit";
                rocheLimitLR.loop = true;
            }
            else
            {
                Debug.LogWarning("Cannot generate Roche limit: no prefab assigned.");
                return;
            }
        }

        if (rocheLimitLR)
        {
            Vector3[] positions = new Vector3[numSamples];
            rocheLimitLR.positionCount = numSamples;
            for (int i = 0; i < numSamples; i++)
            {
                float theta = 2 * Mathf.PI * i / numSamples;
                positions[i] = new Vector3(distance * Mathf.Cos(theta), 0, distance * Mathf.Sin(theta));
            }
            rocheLimitLR.SetPositions(positions);
        }
    }

    public void InstantiatePrefabs(float earthRadius, float lunarRadius, float lunarDistance)
    {
        GenerateEarth(earthRadius);
        GenerateSolidMoon(lunarRadius, lunarDistance);
        DrawRocheLimit(0);

        lights = new Transform[lightPrefabs.Length];
        for (int i = 0; i < lights.Length; i++)
        {
            lights[i] = Instantiate(lightPrefabs[i], transform).transform;
        }
    }

    public void SetLightsVisibility(bool visible)
    {
        foreach (Transform light in lights)
        {
            light.gameObject.SetActive(visible);
        }
    }

    //public void Reset()
    //{
    //    if (moonParticles)
    //    {
    //        // Distribute evenly around a sphere
    //        Vector3 positionCM = Vector3.zero;

    //        float turnFraction = 0.5f * (1 + Mathf.Sqrt(5));  // golden ratio

    //        // Evenly space N bodies around a sphere
    //        for (int i = 0; i < numParticles; i++)
    //        {
    //            Vector3 particlePosition = Vector3.zero;
    //            if (numParticles > 1)
    //            {
    //                float t = i / (numParticles - 1f);
    //                float inclination = Mathf.Acos(1 - 2 * t);
    //                float azimuth = 2 * Mathf.PI * turnFraction * i;

    //                float positionX = Mathf.Sin(inclination) * Mathf.Cos(azimuth);
    //                float positionY = Mathf.Sin(inclination) * Mathf.Sin(azimuth);
    //                float positionZ = Mathf.Cos(inclination);
    //                particlePosition = moonRadius * new Vector3(positionX, positionY, positionZ);
    //            }

    //            moonParticles.GetChild(i).position = particlePosition;
    //            positionCM += particlePosition;
    //        }

    //        positionCM /= numParticles;

    //        // Work in the CM frame (i.e. shift the system to the origin)
    //        for (int i = 0; i < numParticles; i++)
    //        {
    //            moonParticles.GetChild(i).position += moonPosition - positionCM;
    //        }
    //    }

    //    if (phantom)
    //    {
    //        phantom.position = moonPosition;

    //        for (int i = 0; i < numParticles; i++)
    //        {
    //            phantom.GetChild(i).position = moonParticles.GetChild(i).position;
    //        }
    //    }
    //}
}
