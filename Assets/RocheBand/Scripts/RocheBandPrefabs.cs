using UnityEngine;

public class RocheBandPrefabs : MonoBehaviour
{
    [SerializeField] private GameObject primaryPrefab;
    [SerializeField] private GameObject radiusVectorPrefab;
    [SerializeField] private GameObject rigidLimitPrefab;
    [SerializeField] private GameObject fluidLimitPrefab;
    [SerializeField] private GameObject annulusPrefab;
    [SerializeField] private GameObject[] lightPrefabs;

    [HideInInspector] public Transform primary;
    [HideInInspector] public Vector radiusVector;
    [HideInInspector] public LineRenderer rigidLimitLR;
    [HideInInspector] public LineRenderer fluidLimitLR;
    [HideInInspector] public Annulus annulus;
    [HideInInspector] public Transform[] lights;

    private float densityRatio;
    private readonly float rigidFactor = 1.26f;
    private readonly float fluidFactor = 2.44f;

    public void GeneratePrimary(float radius)
    {
        if (primaryPrefab)
        {
            primary = Instantiate(primaryPrefab, transform).transform;
            primary.position = Vector3.zero;
            primary.localScale = 2 * radius * Vector3.one;
            primary.name = "Primary";
        }
        else
        {
            Debug.LogWarning("Cannot generate primary: no prefab assigned.");
        }

        if (radiusVectorPrefab)
        {
            radiusVector = Instantiate(radiusVectorPrefab, transform).GetComponent<Vector>();
            radiusVector.SetPositions(primary.position, primary.position + radius * (Vector3.right + Vector3.up).normalized);
            radiusVector.Redraw();
        }
    }

    public void DrawRigidLimit(float radius, int numSamples = 1000)
    {
        if (!rigidLimitLR)
        {
            if (rigidLimitPrefab)
            {
                rigidLimitLR = Instantiate(rigidLimitPrefab, transform).GetComponent<LineRenderer>();
                rigidLimitLR.transform.position = Vector3.zero;
                rigidLimitLR.name = "Rigid Limit";
                rigidLimitLR.loop = true;
            }
            else
            {
                Debug.LogWarning("Cannot generate rigid limit: no prefab assigned.");
                return;
            }
        }

        if (rigidLimitLR)
        {
            Vector3[] positions = new Vector3[numSamples];
            rigidLimitLR.positionCount = numSamples;
            for (int i = 0; i < numSamples; i++)
            {
                float theta = 2 * Mathf.PI * i / numSamples;
                positions[i] = new Vector3(radius * Mathf.Cos(theta), 0, radius * Mathf.Sin(theta));
            }
            rigidLimitLR.SetPositions(positions);
        }
    }

    public void DrawFluidLimit(float radius, int numSamples = 1000)
    {
        if (!fluidLimitLR)
        {
            if (fluidLimitPrefab)
            {
                fluidLimitLR = Instantiate(fluidLimitPrefab, transform).GetComponent<LineRenderer>();
                fluidLimitLR.transform.position = Vector3.zero;
                fluidLimitLR.name = "Fluid Limit";
                fluidLimitLR.loop = true;
            }
            else
            {
                Debug.LogWarning("Cannot generate fluid limit: no prefab assigned.");
                return;
            }
        }

        if (fluidLimitLR)
        {
            Vector3[] positions = new Vector3[numSamples];
            fluidLimitLR.positionCount = numSamples;
            for (int i = 0; i < numSamples; i++)
            {
                float theta = 2 * Mathf.PI * i / numSamples;
                positions[i] = new Vector3(radius * Mathf.Cos(theta), 0, radius * Mathf.Sin(theta));
            }
            fluidLimitLR.SetPositions(positions);
        }
    }

    public void DrawAnnulus(float innerRadius, float outerRadius)
    {
        if (!annulus)
        {
            if (annulusPrefab)
            {
                annulus = Instantiate(annulusPrefab, transform).GetComponent<Annulus>();
                annulus.transform.position = Vector3.zero;
                annulus.name = "Annulus";
            }
            else
            {
                Debug.LogWarning("Cannot generate annulus: no prefab assigned.");
                return;
            }
        }

        if (annulus)
        {
            annulus.innerRadius = innerRadius;
            annulus.outerRadius = outerRadius;
            annulus.angle = 360;
            annulus.displayPlane = Annulus.DisplayPlane.XZ;
            annulus.Generate();
        }
    }

    public void InstantiateAllPrefabs(float primaryRadius, float densityRatio)
    {
        GeneratePrimary(primaryRadius);

        this.densityRatio = densityRatio;
        float product = Mathf.Pow(this.densityRatio, 0.33f) * primaryRadius;
        DrawRigidLimit(rigidFactor * product);
        DrawFluidLimit(fluidFactor * product);
        DrawAnnulus(rigidFactor * product, fluidFactor * product);

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

    public void SetPrimaryRadius(float value)
    {
        if (primary)
        {
            primary.localScale = 2 * value * Vector3.one;

            if (radiusVector)
            {
                radiusVector.SetPositions(primary.position, primary.position + value * (Vector3.right + Vector3.up).normalized);
                radiusVector.Redraw();
                radiusVector.SetLabelVisibility(value > 1.2f);
            }
        }

        float product = Mathf.Pow(densityRatio, 0.33f) * value;
        DrawRigidLimit(rigidFactor * product);
        DrawFluidLimit(fluidFactor * product);
        DrawAnnulus(rigidFactor * product, fluidFactor * product);
    }
}
