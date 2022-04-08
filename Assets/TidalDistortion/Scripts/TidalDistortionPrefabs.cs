using System.Collections.Generic;
using UnityEngine;

public class TidalDistortionPrefabs : MonoBehaviour
{
    [SerializeField] private GameObject primaryPrefab;
    [SerializeField] private GameObject satellitePrefab;
    [SerializeField] private GameObject rocheLimitPrefab;
    [SerializeField] private GameObject[] lightPrefabs;

    [HideInInspector] public Transform primary;
    [HideInInspector] public Transform satellite;
    [HideInInspector] public Transform particles;
    [HideInInspector] public SphereToEllipsoid satelliteSquasher;
    [HideInInspector] public LineRenderer rocheLimitLR;
    [HideInInspector] public Transform[] lights;

    private float primaryRadius;
    private float primaryDensity;
    private float satelliteDensity;
    private float rocheLimit;

    public void GeneratePrimary(float radius, float density)
    {
        if (primaryPrefab)
        {
            primary = Instantiate(primaryPrefab, transform).transform;
            primary.position = Vector3.zero;
            primary.localScale = 2 * radius * Vector3.one;
            primary.name = "Primary";

            primaryRadius = radius;
            primaryDensity = density;
        }
        else
        {
            Debug.LogWarning("Cannot generate primary: no prefab assigned.");
        }
    }

    public void GenerateSatellite(float distance, float radius, float density)
    {
        if (satellitePrefab)
        {
            satellite = Instantiate(satellitePrefab, transform).transform;
            satellite.position = distance * Vector3.right;
            satellite.localScale = 2 * radius * Vector3.one;
            satellite.name = "Satellite";

            satelliteDensity = density;

            if (!satellite.TryGetComponent(out satelliteSquasher))
            {
                Debug.LogWarning("No SphereToEllipsoid component found on satellite.");
            }
        }
        else
        {
            Debug.LogWarning("Cannot generate satellite: no prefab assigned.");
        }
    }

    public void GenerateParticles()
    {
        if (satellitePrefab && satellite)
        {
            particles = new GameObject("Particles").transform;
            particles.SetParent(transform);

            Mesh mesh = satellite.GetComponent<MeshFilter>().mesh;
            List<Vector3> vertices = new List<Vector3>();
            mesh.GetVertices(vertices);
            for (int i = 0; i < vertices.Count; i++)
            {
                Transform particle = Instantiate(satellitePrefab, particles).transform;
                particle.localScale = 0.12f * Vector3.one;
                particle.position = satellite.position + 0.9f * vertices[i].magnitude * vertices[i].normalized;
            }

            particles.gameObject.SetActive(false);
        }
    }

    public void DrawRocheLimit(int numSamples = 1000)
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
            // Compute Roche limit
            float densityRatio = primaryDensity / satelliteDensity;
            rocheLimit = 2.44f * Mathf.Pow(densityRatio, 0.33f) * primaryRadius;

            Vector3[] positions = new Vector3[numSamples];
            rocheLimitLR.positionCount = numSamples;
            for (int i = 0; i < numSamples; i++)
            {
                float theta = 2 * Mathf.PI * i / numSamples;
                positions[i] = new Vector3(rocheLimit * Mathf.Cos(theta), 0, rocheLimit * Mathf.Sin(theta));
            }
            rocheLimitLR.SetPositions(positions);
        }
    }

    public void InstantiateAllPrefabs(float primaryRadius, float primaryDensity,
        float satelliteDistance, float satelliteRadius, float satelliteDensity)
    {
        GeneratePrimary(primaryRadius, primaryDensity);
        GenerateSatellite(satelliteDistance, satelliteRadius, satelliteDensity);
        GenerateParticles();
        DrawRocheLimit();
        DistortSatellite();

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

    public void SetPrimaryRadius(float value)
    {
        if (primary)
        {
            primary.localScale = 2 * value * Vector3.one;
        }

        primaryRadius = value;
        DrawRocheLimit();
        DistortSatellite();
    }

    public void SetPrimaryDensity(float density)
    {
        primaryDensity = density;
        DrawRocheLimit();
        DistortSatellite();
    }

    public void SetSatelliteDensity(float density)
    {
        satelliteDensity = density;
        DrawRocheLimit();
        DistortSatellite();
    }

    public void DistortSatellite()
    {
        if (satelliteSquasher)
        {
            satelliteSquasher.ShearXY(0.01f * rocheLimit, 0);

            Mesh mesh = satellite.GetComponent<MeshFilter>().mesh;
            List<Vector3> vertices = new List<Vector3>();
            mesh.GetVertices(vertices);
            for (int i = 0; i < vertices.Count; i++)
            {
                particles.GetChild(i).position = satellite.position + 0.9f * vertices[i].magnitude * vertices[i].normalized;
            }

            MeshRenderer renderer = satellite.GetComponent<MeshRenderer>();
            float distance = satellite.position.magnitude;

            if (rocheLimit >= distance && renderer.enabled)
            {
                renderer.enabled = false;
                particles.gameObject.SetActive(true);
            }

            if (rocheLimit < distance && !renderer.enabled)
            {
                renderer.enabled = true;
                particles.gameObject.SetActive(false);
            }
        }
    }
}
