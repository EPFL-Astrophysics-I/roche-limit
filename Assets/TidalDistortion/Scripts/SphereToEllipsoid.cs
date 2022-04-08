using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshFilter))]
public class SphereToEllipsoid : MonoBehaviour
{
    [Range(0, 1)] public float ellipticityXZ = 0;
    [Range(0, 360)] public float angleXZ = 0;
    [Range(0, 1)] public float ellipticityXY = 0;
    [Range(0, 360)] public float angleXY = 0;

    private Mesh mesh;
    private List<Vector3> vertices;
    private List<Vector3> newVertices;

    private void Awake()
    {
        // Get the sphere vertices
        mesh = GetComponent<MeshFilter>().mesh;
        vertices = new List<Vector3>();
        mesh.GetVertices(vertices);

        // Allocate space for deformed sphere vertices
        newVertices = new List<Vector3>(vertices.Count);
        foreach (Vector3 vertex in vertices)
        {
            newVertices.Add(vertex);
        }
    }

    //private void Start()
    //{
    //    UpdateMesh();
    //}

    //private void OnValidate()
    //{
    //    UpdateMesh();
    //}

    //private void UpdateMesh()
    //{
    //    if (mesh == null || vertices == null)
    //    {
    //        return;
    //    }

    //    float phi = angleXZ * Mathf.Deg2Rad;
    //    float theta = angleXY * Mathf.Deg2Rad;
    //    Vector4 col1 = new Vector4(1, 0, 0, 0);
    //    col1 += new Vector4(ellipticityXZ * Mathf.Cos(phi), 0, ellipticityXZ * Mathf.Sin(phi), 0);
    //    col1 += new Vector4(ellipticityXY * Mathf.Cos(theta), ellipticityXY * Mathf.Sin(theta), 0, 0);
    //    Vector4 col2 = new Vector4(0, 1, 0, 0);
    //    col2 += new Vector4(ellipticityXY * Mathf.Sin(theta), -ellipticityXY * Mathf.Cos(theta), 0, 0);
    //    Vector4 col3 = new Vector4(0, 0, 1, 0);
    //    col3 += new Vector4(ellipticityXZ * Mathf.Sin(phi), 0, -ellipticityXZ * Mathf.Cos(phi), 0f);
    //    Vector4 col4 = new Vector4(0, 0, 0, 1f);
    //    Matrix4x4 mat = new Matrix4x4(col1, col2, col3, col4);

    //    for (int i = 0; i < vertices.Count; i++)
    //    {
    //        newVertices[i] = mat.MultiplyPoint3x4(vertices[i]);
    //    }
    //    mesh.SetVertices(newVertices);
    //}

    public void ShearXZ(float ellipticity, float angle)
    {
        if (mesh == null || vertices == null)
        {
            return;
        }

        Vector4 col1 = new Vector4(1, 0, 0, 0);
        col1 += new Vector4(ellipticity * Mathf.Cos(angle), 0, ellipticity * Mathf.Sin(angle), 0);
        Vector4 col2 = new Vector4(0, 1, 0, 0);
        Vector4 col3 = new Vector4(0, 0, 1, 0);
        col3 += new Vector4(ellipticity * Mathf.Sin(angle), 0, -ellipticity * Mathf.Cos(angle), 0);
        Vector4 col4 = new Vector4(0, 0, 0, 1);
        Matrix4x4 matrix = new Matrix4x4(col1, col2, col3, col4);

        for (int i = 0; i < vertices.Count; i++)
        {
            newVertices[i] = matrix.MultiplyPoint3x4(vertices[i]);
        }
        mesh.SetVertices(newVertices);
    }

    public void ShearXY(float ellipticity, float angle)
    {
        if (mesh == null || vertices == null)
        {
            return;
        }

        Vector4 col1 = new Vector4(1, 0, 0, 0);
        col1 += new Vector4(ellipticity * Mathf.Cos(angle), ellipticity * Mathf.Sin(angle), 0, 0);
        Vector4 col2 = new Vector4(0, 1, 0, 0);
        col2 += new Vector4(ellipticity * Mathf.Sin(angle), -ellipticity * Mathf.Cos(angle), 0, 0);
        Vector4 col3 = new Vector4(0, 0, 1, 0);
        Vector4 col4 = new Vector4(0, 0, 0, 1);
        Matrix4x4 matrix = new Matrix4x4(col1, col2, col3, col4);

        for (int i = 0; i < vertices.Count; i++)
        {
            newVertices[i] = matrix.MultiplyPoint3x4(vertices[i]);
        }
        mesh.SetVertices(newVertices);
    }
}
