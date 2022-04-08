using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class Annulus : MonoBehaviour
{
    public bool autoUpdate;

    public enum DisplayPlane { XY, YZ, XZ }
    public DisplayPlane displayPlane = DisplayPlane.XY;

    public float innerRadius = 1;
    public float outerRadius = 3;
    [Range(0, 360)] public float angle = 0;

    private Mesh mesh;
    private List<Vector3> vertexBuffer;
    private List<int> triangleBuffer;
    private int currentVertexIndex;

    public int PositionCount => currentVertexIndex;

    public void Clear()
    {
        if (mesh)
        {
            mesh.Clear();
            mesh = null;
        }

        if (vertexBuffer != null)
        {
            vertexBuffer.Clear();
            vertexBuffer = null;
        }

        if (triangleBuffer != null)
        {
            triangleBuffer.Clear();
            triangleBuffer = null;
        }

        currentVertexIndex = 0;
    }

    public void Generate()
    {
        Clear();

        if (angle < 1)
        {
            return;
        }

        GetComponent<MeshFilter>().mesh = mesh = new Mesh();
        mesh.name = "Annulus Mesh";

        Vector3 iHat = Vector3.right;
        Vector3 jHat = Vector3.up;
        if (displayPlane == DisplayPlane.YZ)
        {
            iHat = Vector3.up;
            jHat = Vector3.forward;
        }
        if (displayPlane == DisplayPlane.XZ)
        {
            jHat = Vector3.forward;
        }

        int numRays = Mathf.Min(360, Mathf.FloorToInt(angle) + 1);
        int numVertices = 2 * numRays;
        vertexBuffer = new List<Vector3>(numVertices);
        int numTriangles = 2 * (numRays - 1);
        if (angle == 360)
        {
            numTriangles += 2;
        }
        triangleBuffer = new List<int>(3 * numTriangles);

        // Always add ray at theta = 0
        Vector3 point1 = innerRadius * iHat;
        Vector3 point2 = outerRadius * iHat;
        vertexBuffer.Add(point1);
        vertexBuffer.Add(point2);

        // Add other points and triangles
        for (int i = 1; i < numRays; i++)
        {
            float theta = i * Mathf.PI / 180f;
            point1 = innerRadius * (Mathf.Cos(theta) * iHat + Mathf.Sin(theta) * jHat);
            point2 = outerRadius * (Mathf.Cos(theta) * iHat + Mathf.Sin(theta) * jHat);
            vertexBuffer.Add(point1);
            vertexBuffer.Add(point2);

            triangleBuffer.Add(i * 2 - 2);
            triangleBuffer.Add(i * 2 + 1);
            triangleBuffer.Add(i * 2 - 1);
            triangleBuffer.Add(i * 2 - 2);
            triangleBuffer.Add(i * 2);
            triangleBuffer.Add(i * 2 + 1);
        }

        if (angle == 360)
        {
            triangleBuffer.Add(numVertices - 2);
            triangleBuffer.Add(1);
            triangleBuffer.Add(numVertices - 1);
            triangleBuffer.Add(numVertices - 2);
            triangleBuffer.Add(0);
            triangleBuffer.Add(1);
        }

        mesh.Clear();
        mesh.SetVertices(vertexBuffer);
        mesh.SetTriangles(triangleBuffer, 0);

        mesh.RecalculateNormals();
        mesh.Optimize();
    }
}
