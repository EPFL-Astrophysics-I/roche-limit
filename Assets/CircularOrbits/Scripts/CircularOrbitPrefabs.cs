using UnityEngine;

public class CircularOrbitPrefabs : MonoBehaviour
{
    [SerializeField] private GameObject primaryPrefab;
    [SerializeField] private GameObject particlePrefab;
    [SerializeField] private GameObject satellitePrefab;
    [SerializeField] private GameObject rocheLimitPrefab;
    [SerializeField] private GameObject rocheVectorPrefab;
    [SerializeField] private GameObject[] lightPrefabs;

    [HideInInspector] public Transform ring;
    [HideInInspector] public Transform satellite;
    [HideInInspector] public LineRenderer rocheLimitLR;
    [HideInInspector] public Vector rocheVector;
    [HideInInspector] public Transform[] lights;


    public void GeneratePrimary(float radius)
    {
        if (primaryPrefab)
        {
            Transform primary = Instantiate(primaryPrefab, transform).transform;
            primary.position = Vector3.zero;
            primary.localScale = 2 * radius * Vector3.one;
            primary.name = "Primary";
        }
        else
        {
            Debug.LogWarning("Cannot generate primary: no prefab assigned.");
        }
    }

    public void GenerateRing(int numParticles, float ringRadius, float ringDispersion, float particleRadius)
    {
        if (particlePrefab)
        {
            ring = new GameObject("Ring").transform;
            ring.SetParent(transform);

            Random.InitState(42);
            float turnFraction = 0.5f * (1 + Mathf.Sqrt(5));  // golden ratio

            for (int i = 0; i < numParticles; i++)
            {
                // Create a particle
                Transform particle = Instantiate(particlePrefab, ring).transform;
                particle.name = "Particle " + i;

                // Space N bodies around a circle
                float radius = Mathf.Clamp(Utils.Random.NormalValue(ringRadius, ringDispersion), ringRadius - 1.5f, ringRadius + 1.5f);
                float theta = 2 * Mathf.PI * turnFraction * i;
                float positionX = radius * Mathf.Cos(theta);
                float positionY = 0;
                float positionZ = radius * Mathf.Sin(theta);

                particle.position = new Vector3(positionX, positionY, positionZ);
                particle.localScale = 2 * particleRadius * Vector3.one;
            }
        }
        else
        {
            Debug.LogWarning("Cannot generate ring: no particle prefab assigned.");
        }
    }

    public void GenerateSatellite(float distance, float radius)
    {
        if (satellitePrefab)
        {
            satellite = Instantiate(satellitePrefab, transform).transform;
            satellite.position = distance * Vector3.right;
            satellite.localScale = 2 * radius * Vector3.one;
            satellite.name = "Satellite";
        }
        else
        {
            Debug.LogWarning("Cannot generate satellite: no prefab assigned.");
        }
    }

    public void GenerateRocheLimit(float radius)
    {
        if (rocheLimitPrefab)
        {
            rocheLimitLR = Instantiate(rocheLimitPrefab, transform).GetComponent<LineRenderer>();
            rocheLimitLR.transform.position = Vector3.zero;
            rocheLimitLR.name = "Roche Limit";

            int numSamples = 1000;
            Vector3[] positions = new Vector3[numSamples];
            rocheLimitLR.positionCount = numSamples;
            for (int i = 0; i < numSamples; i++)
            {
                float theta = 2 * Mathf.PI * i / numSamples;
                positions[i] = new Vector3(radius * Mathf.Cos(theta), 0, radius * Mathf.Sin(theta));
            }
            rocheLimitLR.SetPositions(positions);
            rocheLimitLR.loop = true;
        }
        else
        {
            Debug.LogWarning("Cannot generate Roche limit: no prefab assigned.");
        }

        if (rocheVectorPrefab)
        {
            rocheVector = Instantiate(rocheVectorPrefab, transform).GetComponent<Vector>();
            rocheVector.transform.position = Vector3.zero;
            rocheVector.name = "Roche Vector";

            rocheVector.SetPositions(Vector3.zero, radius * Vector3.right);
            rocheVector.Redraw();
        }
        else
        {
            Debug.LogWarning("Cannot generate Roche vector: no prefab assigned.");
        }
    }

    public void InstantiateAllPrefabs(int numParticles, float primaryRadius, float ringRadius, float ringDispersion,
        float particleRadius, float satelliteDistance, float satelliteRadius, float rocheRadius)
    {
        GeneratePrimary(primaryRadius);
        GenerateRing(numParticles, ringRadius, ringDispersion, particleRadius);
        GenerateSatellite(satelliteDistance, satelliteRadius);
        GenerateRocheLimit(rocheRadius);

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

    public void SetRocheLimitVisibility(bool visible)
    {
        if (rocheLimitLR)
        {
            rocheLimitLR.gameObject.SetActive(visible);
        }
    }
}
