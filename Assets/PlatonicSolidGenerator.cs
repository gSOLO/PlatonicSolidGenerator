#if UNITY_EDITOR
using UnityEditor;  // Needed for drawing labels with Handles in the Unity Editor
#endif
using UnityEngine;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using System;

/// <summary>
/// The available Platonic solids (or 2D pseudo‑solids) to display.
/// </summary>
public enum PlatonicSolid
{
    Tetrahedron,  // A solid with 4 triangular faces
    Cube,         // A solid with 6 square faces
    Octahedron,   // A solid with 8 triangular faces
    Icosahedron,  // A solid with 20 triangular faces
    Dodecahedron, // A solid with 12 pentagonal faces
    Square,       // A 2D shape: a single square face (from a cube)
    Circle        // A 2D shape: a circle defined by evenly spaced vertices
}

/// <summary>
/// Combined connection ordering enum. The first five values use legacy coordinate‑ordering;
/// the remaining four produce a spiral ordering using polar coordinates and a weighted proximity term.
/// </summary>
public enum ConnectionSetting
{
    TopToBottom,
    BottomToTop,
    LeftToRight,
    RightToLeft,
    Random,
    ClockwiseOut,
    ClockwiseIn,
    CounterclockwiseOut,
    CounterclockwiseIn
}

#region PlatonicSolidGizmo Class

/// <summary>
/// PlatonicSolidGizmo is a MonoBehaviour that renders a Platonic (or pseudo‑Platonic) solid in the Unity Editor.
/// It draws a primary wireframe and a secondary wireframe, where the secondary wireframe
/// is rendered as a filled mesh with independently specified face fill and edge colors.
/// Additionally, it generates a cloud of sphere points on the solid's faces, maps letters to these points,
/// processes an input string (activeLettersFilter) via various obfuscation methods, and draws connecting cylinders
/// between the "active" sphere points. Labels can be drawn at each sphere and below the solid.
/// </summary>
public class PlatonicSolidGizmo : MonoBehaviour
{
    #region Public Variables

    // ----- Solid Settings -----
    [Header("Solid Settings")]
    [Tooltip("Choose which Platonic solid to display")]
    public PlatonicSolid solid = PlatonicSolid.Tetrahedron;
    [Tooltip("Uniform scale for the solid")]
    public float scale = 1f;  // Scale applied to the primary solid's vertices

    // ----- Primary Wireframe Options -----
    [Header("Wireframe Options")]
    [Tooltip("Toggle primary wireframe drawing on or off")]
    public bool drawWireframe = true;  // Enables/disables drawing the primary wireframe
    [Tooltip("Color of the primary wireframe of the Platonic solid")]
    public Color wireframeColor = Color.green;

    // ----- Secondary Wireframe Options -----
    [Header("Secondary Wireframe Options")]
    [Tooltip("Toggle drawing of the secondary (solid) wireframe")]
    public bool drawSecondaryWireframe = true;  // Enables secondary wireframe drawing
    [Tooltip("Uniform scale for the secondary wireframe (independent of the main solid scale)")]
    public float secondaryScale = 1.0f;  // Scale factor for the secondary wireframe
    [Tooltip("Fill color for the secondary solid faces")]
    public Color secondaryFaceFillColor = new Color(1f, 1f, 1f, 0.5f);  // Face fill color
    [Tooltip("Color for the secondary wireframe edges")]
    public Color secondaryWireframeEdgeColor = Color.black;  // Edge (line) color for secondary wireframe

    // ----- Sphere Cloud Settings -----
    [Header("Sphere Cloud Settings")]
    [Tooltip("Exact total number of spheres to distribute over the faces")]
    public int spherePointCount = 50;  // Number of sphere points to generate
    [Tooltip("Base radius for the first sphere")]
    public float sphereRadius = 0.05f;
    [Tooltip("Additional radius added per sphere (each subsequent sphere is larger)")]
    public float sphereSizeIncrement = 0.01f;

    // ----- Sphere Color Gradient -----
    [Header("Sphere Color Gradient")]
    [Tooltip("Color (with opacity) for the first sphere")]
    public Color startSphereColor = new Color(1f, 0f, 0f, 1f);
    [Tooltip("Color (with opacity) for the last sphere")]
    public Color endSphereColor = new Color(0f, 0f, 1f, 0.5f);
    // A gradient is applied along the chain of sphere points.

    // ----- Letter Mapping Settings -----
    [Header("Letter Mapping Settings")]
    [Tooltip("If false, letters are assigned in A–Z order; if true, use alternate ordering")]
    public bool useAlternateAlphabet = false;
    [Tooltip("Alphabet used for A–Z mapping")]
    public string alphabetAZ = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    [Tooltip("Alternate alphabet ordering")]
    public string alphabetEQ = "ESIARNTOLCDUGPMHBYFVKWZXJQ";
    [Tooltip("Active letter chain. Type a sentence here which will be processed (via obfuscation methods) " +
             "to produce the chain of active letters. For example, 'AGPA' selects spheres with A, then G, then P, then A (even if letters repeat)")]
    public string activeLettersFilter = "";  // User input for active chain filtering

    // ----- Active Letters Obfuscation Settings -----
    [Header("Active Letters Obfuscation Settings")]
    // Toggles for string transformation methods.
    [Tooltip("Toggle the 'Flip' method (reverses the word order) on or off")]
    public bool obfFlip = false;
    [Tooltip("Toggle the 'Rotate' method (rotates the words) on or off")]
    public bool obfRotate = false;
    [Tooltip("Offset for the Rotate method (optional)")]
    public int obfRotateOffset = 0;
    [Tooltip("Toggle the 'CutUp' method (randomly shuffles the words) on or off")]
    public bool obfCutUp = false;
    [Tooltip("Toggle the 'RemoveWhitespace' method on or off")]
    public bool obfRemoveWhitespace = false;
    [Tooltip("Toggle the 'Reverse' method on or off")]
    public bool obfReverse = false;
    [Tooltip("Toggle the 'Shift' method on or off")]
    public bool obfShift = false;
    [Tooltip("Offset for the Shift method (optional)")]
    public int obfShiftOffset = 0;
    [Tooltip("Toggle the 'Disemvowel' method on or off. When enabled, vowels are removed from both the active filter and the alphabets")]
    public bool obfDisemvowel = false;
    [Tooltip("Toggle the 'Squeeze' method on or off")]
    public bool obfSqueeze = false;
    [Tooltip("Toggle the 'Deduplicate' method on or off")]
    public bool obfDeduplicate = false;
    [Tooltip("Toggle the 'Rearrange' method on or off")]
    public bool obfRearrange = false;

    // ----- Endpoint Options -----
    [Header("Endpoint Options")]
    [Tooltip("When enabled, only the first and last spheres are drawn (the last is drawn as a cube). " +
             "Intermediate sphere points are not drawn, but labels are still placed at these points.")]
    public bool drawEndpointsOnly = false;

    // ----- Label Options -----
    [Header("Label Options")]
    [Tooltip("Toggle all labels (sphere labels, processed filter label, etc.) on or off")]
    public bool drawLabels = true;

    // ----- Candidate Generation Settings -----
    [Header("Candidate Generation Settings")]
    [Tooltip("Subdivision level for candidate point generation on each face (for Circle, clamped to determine vertex count)")]
    public int candidateSubdivision = 10;

    // ----- Cylinder Settings -----
    [Header("Cylinder Settings")]
    [Tooltip("Radius for the connecting cylinders between sphere centers")]
    public float cylinderRadius = 0.05f;

    // ----- Connection Ordering Settings -----
    [Header("Connection Ordering Settings")]
    [Tooltip("Combined connection setting that determines how sphere points are ordered. " +
             "Legacy options (TopToBottom, etc.) sort by coordinate; spiral options (ClockwiseOut/In, CounterclockwiseOut/In) produce a spiral. " +
             "The spiral ordering uses both angle and proximity. Adjust the proximity influence with the variable below.")]
    public ConnectionSetting connectionSetting = ConnectionSetting.TopToBottom;
    [Tooltip("Used only for legacy ordering options; ignored for spiral ordering")]
    public int startingVertexIndex = 0;
    [Tooltip("Weight factor (typically between 0 and 1) for the proximity term in spiral ordering. Lower values emphasize angle more.")]
    public float spiralProximityWeight = 0.3f;

    #endregion

    #region Private Variables

    // Cached cylinder mesh for drawing connecting cylinders.
    private static Mesh cylinderMesh;

    #endregion

    #region OnDrawGizmos

    /// <summary>
    /// OnDrawGizmos is called by Unity to render custom gizmos in the Scene view.
    /// It performs the following steps:
    /// 1. Retrieves and transforms primary vertices.
    /// 2. Draws the primary wireframe.
    /// 3. Draws the secondary wireframe (first filling faces, then drawing edge lines).
    /// 4. Generates candidate sphere points on the solid’s faces.
    /// 5. Orders the sphere points based on the connection setting.
    /// 6. Maps letters to each sphere point.
    /// 7. Processes the active letter filter to build an active chain of sphere points.
    /// 8. Assigns a color gradient to the sphere points.
    /// 9. Draws the sphere points (or endpoints only) and labels.
    /// 10. Draws connecting cylinders between consecutive visible sphere points.
    /// 11. Draws a label below the solid displaying the processed active filter.
    /// </summary>
    private void OnDrawGizmos()
    {
        // --- Step 1: Retrieve and transform primary vertices ---
        Vector3[] vertices = GetVertices(solid);
        if (vertices == null || vertices.Length == 0)
            return;
        // Apply primary scale and position.
        for (int i = 0; i < vertices.Length; i++)
            vertices[i] = transform.position + vertices[i] * scale;

        // --- Step 2: Draw the primary wireframe ---
        if (drawWireframe)
        {
            // Determine the smallest edge length for drawing a simple wireframe.
            float minEdgeLength = float.MaxValue;
            for (int i = 0; i < vertices.Length; i++)
            {
                for (int j = i + 1; j < vertices.Length; j++)
                {
                    float d = Vector3.Distance(vertices[i], vertices[j]);
                    if (d > 0 && d < minEdgeLength)
                        minEdgeLength = d;
                }
            }
            float tol = 0.01f;
            Gizmos.color = wireframeColor;
            // Draw lines between vertices that are approximately the minimum edge length apart.
            for (int i = 0; i < vertices.Length; i++)
            {
                for (int j = i + 1; j < vertices.Length; j++)
                {
                    float d = Vector3.Distance(vertices[i], vertices[j]);
                    if (Mathf.Abs(d - minEdgeLength) < tol)
                        Gizmos.DrawLine(vertices[i], vertices[j]);
                }
            }
        }

        // --- Step 3: Draw the secondary wireframe (filled faces + edge lines) ---
        if (drawSecondaryWireframe)
        {
            // Retrieve secondary vertices with independent scaling.
            Vector3[] secondaryVertices = GetVertices(solid);
            if (secondaryVertices != null && secondaryVertices.Length > 0)
            {
                // Apply secondary scale and position.
                for (int i = 0; i < secondaryVertices.Length; i++)
                    secondaryVertices[i] = transform.position + secondaryVertices[i] * secondaryScale;

                // Get the face definitions for the secondary wireframe.
                List<int[]> solidFaces = GetFaces(solid);

                // Generate a solid mesh (for face fill) from the secondary vertices.
                Mesh solidMesh = GenerateSolidMesh(secondaryVertices, solidFaces);
                if (solidMesh != null)
                {
                    // For Square and Circle, make the mesh double-sided.
                    if (solid == PlatonicSolid.Square || solid == PlatonicSolid.Circle)
                    {
                        solidMesh = MakeMeshDoubleSided(solidMesh);
                    }
                    // Draw the solid faces using the secondary face fill color.
                    Gizmos.color = secondaryFaceFillColor;
                    Gizmos.DrawMesh(solidMesh);
                }
                // Draw the wireframe edges over the secondary mesh.
                Gizmos.color = secondaryWireframeEdgeColor;
                foreach (var face in solidFaces)
                {
                    // Draw a line loop for each face.
                    for (int i = 0; i < face.Length; i++)
                    {
                        Vector3 v0 = secondaryVertices[face[i]];
                        Vector3 v1 = secondaryVertices[face[(i + 1) % face.Length]];
                        Gizmos.DrawLine(v0, v1);
                    }
                }
            }
        }

        // --- Step 4: Generate candidate sphere points ---
        List<Vector3> finalPoints = new List<Vector3>();
        if (spherePointCount == vertices.Length)
        {
            // If the number of sphere points equals the number of vertices, use them directly.
            finalPoints.AddRange(vertices);
        }
        else
        {
            // Otherwise, generate candidate points based on face subdivision.
            List<int[]> faces = GetFaces(solid);
            if (faces == null || faces.Count == 0)
                return;
            List<Vector3> candidatePoints = new List<Vector3>();
            // For each face, generate candidate points using barycentric subdivision.
            foreach (var faceIndices in faces)
            {
                Vector3[] faceVerts = new Vector3[faceIndices.Length];
                for (int i = 0; i < faceIndices.Length; i++)
                    faceVerts[i] = vertices[faceIndices[i]];
                List<Vector3> faceCandidates = GenerateCandidatePointsOnPolygon(faceVerts, candidateSubdivision);
                candidatePoints.AddRange(faceCandidates);
            }
            candidatePoints = RemoveDuplicatePoints(candidatePoints, 0.0001f);
            if (candidatePoints.Count < spherePointCount)
            {
                Debug.LogWarning("Not enough candidate points generated. Increase candidateSubdivision.");
                spherePointCount = candidatePoints.Count;
            }
            if (spherePointCount > vertices.Length)
                finalPoints = FarthestPointSamplingWithPreselected(candidatePoints, spherePointCount, new List<Vector3>(vertices));
            else
                finalPoints = FarthestPointSampling(candidatePoints, spherePointCount, vertices);
        }

        // --- Step 5: Order the sphere points for connection ---
        List<Vector3> orderedPoints = new List<Vector3>();
        if (connectionSetting == ConnectionSetting.TopToBottom ||
            connectionSetting == ConnectionSetting.BottomToTop ||
            connectionSetting == ConnectionSetting.LeftToRight ||
            connectionSetting == ConnectionSetting.RightToLeft ||
            connectionSetting == ConnectionSetting.Random)
        {
            // Legacy coordinate-based ordering.
            orderedPoints = new List<Vector3>(finalPoints);
            const float orderTol = 0.001f;
            switch (connectionSetting)
            {
                case ConnectionSetting.TopToBottom:
                    orderedPoints.Sort((a, b) =>
                        (Mathf.Abs(a.y - b.y) < orderTol) ? a.x.CompareTo(b.x) : b.y.CompareTo(a.y));
                    break;
                case ConnectionSetting.BottomToTop:
                    orderedPoints.Sort((a, b) =>
                        (Mathf.Abs(a.y - b.y) < orderTol) ? a.x.CompareTo(b.x) : a.y.CompareTo(b.y));
                    break;
                case ConnectionSetting.LeftToRight:
                    orderedPoints.Sort((a, b) =>
                        (Mathf.Abs(a.x - b.x) < orderTol) ? a.y.CompareTo(b.y) : a.x.CompareTo(b.x));
                    break;
                case ConnectionSetting.RightToLeft:
                    orderedPoints.Sort((a, b) =>
                        (Mathf.Abs(a.x - b.x) < orderTol) ? a.y.CompareTo(b.y) : b.x.CompareTo(a.x));
                    break;
                case ConnectionSetting.Random:
                    // Random shuffle.
                    for (int i = 0; i < orderedPoints.Count; i++)
                    {
                        int rnd = UnityEngine.Random.Range(i, orderedPoints.Count);
                        Vector3 temp = orderedPoints[i];
                        orderedPoints[i] = orderedPoints[rnd];
                        orderedPoints[rnd] = temp;
                    }
                    break;
            }
            // Rotate the ordered list based on startingVertexIndex.
            int startIndex = Mathf.Clamp(startingVertexIndex, 0, orderedPoints.Count - 1);
            List<Vector3> rotatedPoints = new List<Vector3>();
            for (int i = startIndex; i < orderedPoints.Count; i++)
                rotatedPoints.Add(orderedPoints[i]);
            for (int i = 0; i < startIndex; i++)
                rotatedPoints.Add(orderedPoints[i]);
            orderedPoints = rotatedPoints;
        }
        else
        {
            // Spiral ordering using polar coordinates.
            Vector3 center = Vector3.zero;
            foreach (Vector3 pt in finalPoints)
                center += pt;
            center /= finalPoints.Count;

            float maxR = 0f;
            List<(Vector3 pt, float r, float theta)> polarList = new List<(Vector3, float, float)>();
            foreach (Vector3 pt in finalPoints)
            {
                Vector3 d = pt - center;
                float r = d.magnitude;
                if (r > maxR) maxR = r;
                float theta = Mathf.Atan2(d.y, d.x);
                if (theta < 0) theta += 2 * Mathf.PI;
                polarList.Add((pt, r, theta));
            }

            bool clockwise = (connectionSetting == ConnectionSetting.ClockwiseOut || connectionSetting == ConnectionSetting.ClockwiseIn);
            bool spiralOut = (connectionSetting == ConnectionSetting.ClockwiseOut || connectionSetting == ConnectionSetting.CounterclockwiseOut);

            polarList.Sort((a, b) =>
            {
                float angleA = clockwise ? (2 * Mathf.PI - a.theta) : a.theta;
                float angleB = clockwise ? (2 * Mathf.PI - b.theta) : b.theta;
                float normA = (maxR > 0 ? a.r / maxR : 0f);
                float normB = (maxR > 0 ? b.r / maxR : 0f);
                float radiusTermA = spiralOut ? normA : (1f - normA);
                float radiusTermB = spiralOut ? normB : (1f - normB);
                float keyA = angleA + spiralProximityWeight * radiusTermA;
                float keyB = angleB + spiralProximityWeight * radiusTermB;
                return keyA.CompareTo(keyB);
            });

            List<Vector3> sortedByAngle = polarList.Select(tuple => tuple.pt).ToList();
            if (startingVertexIndex > 0 && startingVertexIndex < sortedByAngle.Count)
            {
                List<Vector3> rotated = new List<Vector3>();
                for (int i = startingVertexIndex; i < sortedByAngle.Count; i++)
                    rotated.Add(sortedByAngle[i]);
                for (int i = 0; i < startingVertexIndex; i++)
                    rotated.Add(sortedByAngle[i]);
                sortedByAngle = rotated;
            }
            orderedPoints = sortedByAngle;
        }

        // --- Step 6: Map letters to each ordered sphere point ---
        int pointCount = orderedPoints.Count;
        string[] letterMapping = new string[pointCount];
        for (int i = 0; i < pointCount; i++)
            letterMapping[i] = "";
        // Choose which alphabet to use; optionally remove vowels.
        string effectiveAlphabet = useAlternateAlphabet ? alphabetEQ : alphabetAZ;
        if (obfDisemvowel)
            effectiveAlphabet = effectiveAlphabet.Disemvowel();
        int effectiveAlphabetLength = effectiveAlphabet.Length;
        if (pointCount < effectiveAlphabetLength)
        {
            for (int i = 0; i < effectiveAlphabetLength; i++)
            {
                int idx = i % pointCount;
                letterMapping[idx] += effectiveAlphabet[i % effectiveAlphabetLength];
            }
        }
        else
        {
            for (int i = 0; i < effectiveAlphabetLength; i++)
            {
                int idx = Mathf.RoundToInt(i * (pointCount - 1) / (float)(effectiveAlphabetLength - 1));
                letterMapping[idx] += effectiveAlphabet[i % effectiveAlphabetLength];
            }
        }

        // --- Step 7: Process the active letter filter ---
        string processedFilter = ProcessActiveLettersFilter(activeLettersFilter);

        // --- Step 8: Build the active chain of sphere point indices ---
        List<int> activeChain = new List<int>();
        if (string.IsNullOrEmpty(processedFilter))
        {
            for (int i = 0; i < pointCount; i++)
                activeChain.Add(i);
        }
        else
        {
            int currentIndex = 0;
            for (int f = 0; f < processedFilter.Length; f++)
            {
                char c = processedFilter[f];
                for (int i = currentIndex; i < currentIndex + pointCount; i++)
                {
                    int idx = i % pointCount;
                    if (letterMapping[idx].IndexOf(c) >= 0)
                    {
                        activeChain.Add(idx);
                        currentIndex = idx + 1;
                        break;
                    }
                }
            }
        }

        // --- Step 9: Build lists of visible points and assign a color gradient ---
        List<Vector3> visiblePoints = new List<Vector3>();
        List<Color> visibleColors = new List<Color>();
        if (string.IsNullOrEmpty(activeLettersFilter))
        {
            Color[] sphereColors = new Color[pointCount];
            for (int i = 0; i < pointCount; i++)
            {
                float t = (pointCount > 1) ? (float)i / (pointCount - 1) : 0f;
                sphereColors[i] = Color.Lerp(startSphereColor, endSphereColor, t);
            }
            foreach (int idx in activeChain)
            {
                visiblePoints.Add(orderedPoints[idx]);
                visibleColors.Add(sphereColors[idx]);
            }
        }
        else
        {
            int visibleCount = activeChain.Count;
            for (int i = 0; i < visibleCount; i++)
            {
                float t = (visibleCount > 1) ? (float)i / (visibleCount - 1) : 0f;
                visibleColors.Add(Color.Lerp(startSphereColor, endSphereColor, t));
            }
            foreach (int idx in activeChain)
            {
                visiblePoints.Add(orderedPoints[idx]);
            }
        }

        // --- Step 10: Draw the sphere points and labels ---
        for (int i = 0; i < visiblePoints.Count; i++)
        {
            float r = sphereRadius + (activeChain[i]) * sphereSizeIncrement;
            if (!drawEndpointsOnly)
            {
                Gizmos.color = visibleColors[i];
                Gizmos.DrawSphere(visiblePoints[i], r);
            }
            else
            {
                if (i == 0)
                {
                    Gizmos.color = visibleColors[i];
                    Gizmos.DrawSphere(visiblePoints[i], r);
                }
                else if (i == visiblePoints.Count - 1)
                {
                    Gizmos.color = visibleColors[i];
                    Gizmos.DrawCube(visiblePoints[i], new Vector3(2 * r, 2 * r, 2 * r));
                }
            }
#if UNITY_EDITOR
            if (drawLabels)
                Handles.Label(visiblePoints[i] + Vector3.up * (r * 1.5f), letterMapping[activeChain[i]]);
#endif
        }

        // --- Step 11: Draw connecting cylinders between consecutive visible sphere points ---
        EnsureCylinderMesh();
        for (int i = 0; i < visiblePoints.Count - 1; i++)
        {
            Vector3 p0 = visiblePoints[i];
            Vector3 p1 = visiblePoints[i + 1];
            Vector3 direction = p1 - p0;
            float distance = direction.magnitude;
            Vector3 mid = (p0 + p1) / 2f;
            Quaternion rotation = Quaternion.FromToRotation(Vector3.up, direction);
            Vector3 cylScale = new Vector3(cylinderRadius / 0.5f, distance / 2f, cylinderRadius / 0.5f);
            Gizmos.color = visibleColors[i + 1];
            Gizmos.DrawMesh(cylinderMesh, mid, rotation, cylScale);
        }

#if UNITY_EDITOR
        // --- Step 12: Draw a centered label below the solid displaying the processed active filter ---
        if (drawLabels)
        {
            Vector3 min = vertices[0];
            Vector3 max = vertices[0];
            foreach (Vector3 v in vertices)
            {
                min = Vector3.Min(min, v);
                max = Vector3.Max(max, v);
            }
            Vector3 centerSolid = (min + max) / 2;
            Vector3 labelPos = new Vector3(centerSolid.x, min.y - sphereRadius * 2, centerSolid.z);
            Handles.Label(labelPos, $"Processed Filter: {processedFilter}");
        }
#endif
    }

    #endregion

    #region String Obfuscation Processing

    /// <summary>
    /// Processes the activeLettersFilter string by applying a series of obfuscation methods.
    /// Methods are applied in a fixed order: ToUpper (always), then optionally:
    /// Flip, Rotate, CutUp, RemoveWhitespace, Reverse, Shift, Disemvowel, Squeeze, Deduplicate, and Rearrange.
    /// </summary>
    /// <param name="filter">The original filter string input by the user.</param>
    /// <returns>The transformed (obfuscated) string.</returns>
    private string ProcessActiveLettersFilter(string filter)
    {
        if (string.IsNullOrEmpty(filter))
            return filter;
        // Always convert to uppercase.
        filter = filter.ToUpper();
        if (obfFlip)
            filter = filter.Flip();
        if (obfRotate)
            filter = filter.Rotate(obfRotateOffset);
        if (obfCutUp)
            filter = filter.CutUp();
        if (obfRemoveWhitespace)
            filter = filter.RemoveWhitespace();
        if (obfReverse)
            filter = filter.Reverse();
        if (obfShift)
            filter = filter.Shift(obfShiftOffset);
        if (obfDisemvowel)
            filter = filter.Disemvowel();
        if (obfSqueeze)
            filter = filter.Squeeze();
        if (obfDeduplicate)
            filter = filter.Deduplicate();
        if (obfRearrange)
            filter = filter.Rearrange();
        return filter;
    }

    #endregion

    #region Sampling and Geometry Helpers

    /// <summary>
    /// Caches a cylinder mesh based on Unity's built-in Cylinder primitive.
    /// Creates a temporary Cylinder object to obtain its mesh, then destroys it.
    /// </summary>
    private static void EnsureCylinderMesh()
    {
        if (cylinderMesh == null)
        {
            GameObject temp = GameObject.CreatePrimitive(PrimitiveType.Cylinder);
            cylinderMesh = temp.GetComponent<MeshFilter>().sharedMesh;
            DestroyImmediate(temp);
        }
    }

    /// <summary>
    /// Generates a solid mesh for the secondary wireframe by fan‑triangulating each face.
    /// This mesh is used to fill the faces with a solid color.
    /// </summary>
    /// <param name="vertices">Array of secondary vertices (transformed with secondaryScale and position).</param>
    /// <param name="faces">List of faces (each face is defined by an array of vertex indices).</param>
    /// <returns>A Mesh representing the secondary solid's faces.</returns>
    private Mesh GenerateSolidMesh(Vector3[] vertices, List<int[]> faces)
    {
        Mesh mesh = new Mesh();
        mesh.name = "SecondaryWireframeMesh";
        mesh.vertices = vertices;

        List<int> triangles = new List<int>();

        // For each face, sort the indices and then perform fan‑triangulation.
        foreach (int[] face in faces)
        {
            if (face.Length < 3)
                continue;

            // Sort the face indices to ensure consistent winding.
            int[] sortedFace = SortFaceIndices(face, vertices);

            // Fan-triangulation: from the first vertex, create triangles with subsequent pairs.
            for (int i = 1; i < sortedFace.Length - 1; i++)
            {
                triangles.Add(sortedFace[0]);
                triangles.Add(sortedFace[i]);
                triangles.Add(sortedFace[i + 1]);
            }
        }

        mesh.triangles = triangles.ToArray();
        mesh.RecalculateNormals();
        return mesh;
    }

    /// <summary>
    /// Sorts the vertex indices of a face into a consistent (e.g., counterclockwise) winding order.
    /// Computes the face centroid and a representative normal, then sorts based on the angle between
    /// a reference direction (from the centroid to the first vertex) and each vertex's direction.
    /// </summary>
    /// <param name="face">Array of vertex indices defining the face.</param>
    /// <param name="vertices">Array of vertices.</param>
    /// <returns>An array of vertex indices sorted into a consistent winding order.</returns>
    private int[] SortFaceIndices(int[] face, Vector3[] vertices)
    {
        List<Vector3> faceVerts = new List<Vector3>();
        foreach (int index in face)
            faceVerts.Add(vertices[index]);

        // Compute the centroid of the face.
        Vector3 centroid = Vector3.zero;
        foreach (Vector3 v in faceVerts)
            centroid += v;
        centroid /= faceVerts.Count;

        // Compute a representative normal from the first three vertices.
        Vector3 normal = Vector3.zero;
        if (faceVerts.Count >= 3)
            normal = Vector3.Cross(faceVerts[1] - faceVerts[0], faceVerts[2] - faceVerts[0]).normalized;

        // Use the direction from the centroid to the first vertex as a reference.
        Vector3 refDir = (faceVerts[0] - centroid).normalized;

        List<(int index, float angle)> indexedAngles = new List<(int, float)>();
        for (int i = 0; i < faceVerts.Count; i++)
        {
            Vector3 dir = (faceVerts[i] - centroid).normalized;
            float angle = Mathf.Atan2(Vector3.Dot(Vector3.Cross(refDir, dir), normal),
                                      Vector3.Dot(refDir, dir));
            indexedAngles.Add((face[i], angle));
        }
        indexedAngles.Sort((a, b) => a.angle.CompareTo(b.angle));
        return indexedAngles.Select(pair => pair.index).ToArray();
    }

    /// <summary>
    /// Standard farthest point sampling.
    /// Selects 'count' points from a list of candidate points by iteratively choosing the candidate farthest
    /// from the already selected points. Starts with a candidate near an original vertex, if possible.
    /// </summary>
    /// <param name="candidates">List of candidate points.</param>
    /// <param name="count">Desired number of points.</param>
    /// <param name="originalVertices">Original vertices of the solid (used for seeding).</param>
    /// <returns>A list of evenly distributed points.</returns>
    private List<Vector3> FarthestPointSampling(List<Vector3> candidates, int count, Vector3[] originalVertices)
    {
        List<Vector3> selected = new List<Vector3>();
        int n = candidates.Count;
        if (n == 0 || count <= 0)
            return selected;
        float tol = 0.0001f;
        int initialIndex = -1;
        for (int i = 0; i < n; i++)
        {
            foreach (Vector3 v in originalVertices)
            {
                if (Vector3.Distance(candidates[i], v) < tol)
                {
                    initialIndex = i;
                    break;
                }
            }
            if (initialIndex != -1)
                break;
        }
        if (initialIndex == -1)
            initialIndex = UnityEngine.Random.Range(0, n);
        selected.Add(candidates[initialIndex]);
        float[] minDist = new float[n];
        for (int i = 0; i < n; i++)
            minDist[i] = Vector3.Distance(candidates[i], selected[0]);
        while (selected.Count < count)
        {
            int bestIndex = -1;
            float bestDist = -1f;
            for (int i = 0; i < n; i++)
            {
                if (selected.Contains(candidates[i]))
                    continue;
                if (minDist[i] > bestDist)
                {
                    bestDist = minDist[i];
                    bestIndex = i;
                }
            }
            if (bestIndex == -1)
                break;
            Vector3 bestCandidate = candidates[bestIndex];
            selected.Add(bestCandidate);
            for (int i = 0; i < n; i++)
            {
                float d = Vector3.Distance(candidates[i], bestCandidate);
                if (d < minDist[i])
                    minDist[i] = d;
            }
        }
        return selected;
    }

    /// <summary>
    /// Modified farthest point sampling that includes a preselected set of points,
    /// then fills in additional points using farthest point sampling.
    /// </summary>
    /// <param name="candidates">List of candidate points.</param>
    /// <param name="count">Total desired number of points.</param>
    /// <param name="preselected">Points that must be included.</param>
    /// <returns>A list containing the preselected points plus additional evenly distributed points.</returns>
    private List<Vector3> FarthestPointSamplingWithPreselected(List<Vector3> candidates, int count, List<Vector3> preselected)
    {
        List<Vector3> selected = new List<Vector3>(preselected);
        int n = candidates.Count;
        if (n == 0 || count <= selected.Count)
            return selected;
        float tol = 0.0001f;
        float[] minDist = new float[n];
        for (int i = 0; i < n; i++)
        {
            float dmin = float.MaxValue;
            foreach (Vector3 s in selected)
            {
                float d = Vector3.Distance(candidates[i], s);
                if (d < dmin)
                    dmin = d;
            }
            minDist[i] = dmin;
        }
        while (selected.Count < count)
        {
            int bestIndex = -1;
            float bestDist = -1f;
            for (int i = 0; i < n; i++)
            {
                bool alreadySelected = false;
                foreach (Vector3 s in selected)
                {
                    if (Vector3.Distance(candidates[i], s) < tol)
                    {
                        alreadySelected = true;
                        break;
                    }
                }
                if (alreadySelected)
                    continue;
                if (minDist[i] > bestDist)
                {
                    bestDist = minDist[i];
                    bestIndex = i;
                }
            }
            if (bestIndex == -1)
                break;
            Vector3 bestCandidate = candidates[bestIndex];
            selected.Add(bestCandidate);
            for (int i = 0; i < n; i++)
            {
                float d = Vector3.Distance(candidates[i], bestCandidate);
                if (d < minDist[i])
                    minDist[i] = d;
            }
        }
        return selected;
    }

    /// <summary>
    /// Generates candidate points on a polygon using fan‑triangulation.
    /// The polygon is subdivided into triangles, and a barycentric grid is generated on each triangle.
    /// </summary>
    /// <param name="polyVerts">Array of vertices defining the polygon.</param>
    /// <param name="subdivision">Number of subdivisions per triangle.</param>
    /// <returns>A list of candidate points over the polygon.</returns>
    private List<Vector3> GenerateCandidatePointsOnPolygon(Vector3[] polyVerts, int subdivision)
    {
        List<Vector3> points = new List<Vector3>();
        if (polyVerts.Length < 3)
            return points;
        int triangleCount = polyVerts.Length - 2;
        for (int i = 0; i < triangleCount; i++)
        {
            Vector3 A = polyVerts[0];
            Vector3 B = polyVerts[i + 1];
            Vector3 C = polyVerts[i + 2];
            for (int j = 0; j <= subdivision; j++)
            {
                for (int k = 0; k <= subdivision - j; k++)
                {
                    float u = j / (float)subdivision;
                    float v = k / (float)subdivision;
                    float w = 1f - u - v;
                    if (w < 0)
                        continue;
                    Vector3 p = A * w + B * u + C * v;
                    points.Add(p);
                }
            }
        }
        return points;
    }

    /// <summary>
    /// Removes duplicate points from a list based on a specified tolerance.
    /// </summary>
    /// <param name="points">List of points.</param>
    /// <param name="tolerance">Distance under which points are considered duplicates.</param>
    /// <returns>A list of unique points.</returns>
    private List<Vector3> RemoveDuplicatePoints(List<Vector3> points, float tolerance)
    {
        List<Vector3> unique = new List<Vector3>();
        foreach (var p in points)
        {
            bool duplicate = false;
            foreach (var q in unique)
            {
                if (Vector3.Distance(p, q) < tolerance)
                {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate)
                unique.Add(p);
        }
        return unique;
    }

    /// <summary>
    /// Returns the vertices for the selected solid.
    /// For 3D solids, vertices are defined symmetrically.
    /// For 2D shapes (Square and Circle), vertices lie in the XY plane.
    /// </summary>
    /// <param name="solidType">The selected Platonic solid.</param>
    /// <returns>An array of vertices for the solid.</returns>
    private Vector3[] GetVertices(PlatonicSolid solidType)
    {
        switch (solidType)
        {
            case PlatonicSolid.Tetrahedron:
                return new Vector3[]
                {
                    new Vector3( 1,  1,  1),
                    new Vector3( 1, -1, -1),
                    new Vector3(-1,  1, -1),
                    new Vector3(-1, -1,  1)
                };
            case PlatonicSolid.Cube:
                return new Vector3[]
                {
                    new Vector3(-1, -1, -1),
                    new Vector3(-1, -1,  1),
                    new Vector3(-1,  1, -1),
                    new Vector3(-1,  1,  1),
                    new Vector3( 1, -1, -1),
                    new Vector3( 1, -1,  1),
                    new Vector3( 1,  1, -1),
                    new Vector3( 1,  1,  1)
                };
            case PlatonicSolid.Octahedron:
                return new Vector3[]
                {
                    new Vector3( 1,  0,  0),
                    new Vector3(-1,  0,  0),
                    new Vector3( 0,  1,  0),
                    new Vector3( 0, -1,  0),
                    new Vector3( 0,  0,  1),
                    new Vector3( 0,  0, -1)
                };
            case PlatonicSolid.Icosahedron:
                {
                    float a = 1.0f;
                    float b = 0.0f;
                    float c = (1f + Mathf.Sqrt(5f)) / 2f; // phi
                    return new Vector3[]
                    {
                        new Vector3(+b, +a, +c),
                        new Vector3(+b, +a, -c),
                        new Vector3(+b, -a, +c),
                        new Vector3(+b, -a, -c),

                        new Vector3(+a, +c, +b),
                        new Vector3(+a, -c, +b),
                        new Vector3(-a, +c, +b),
                        new Vector3(-a, -c, +b),

                        new Vector3(+c, +b, +a),
                        new Vector3(+c, +b, -a),
                        new Vector3(-c, +b, +a),
                        new Vector3(-c, +b, -a)
                    };
                }
            case PlatonicSolid.Dodecahedron:
                {
                    float a = 1.0f;
                    float b = 0.0f;
                    float c = (1f + Mathf.Sqrt(5f)) / 2f; // phi
                    float d = 1f / c; // invPhi
                    List<Vector3> verts = new List<Vector3>();
                    verts.Add(new Vector3(+a, +a, +a));
                    verts.Add(new Vector3(+a, +a, -a));
                    verts.Add(new Vector3(+a, -a, +a));
                    verts.Add(new Vector3(-a, +a, +a));

                    verts.Add(new Vector3(+a, -a, -a));
                    verts.Add(new Vector3(-a, +a, -a));
                    verts.Add(new Vector3(-a, -a, +a));
                    verts.Add(new Vector3(-a, -a, -a));

                    verts.Add(new Vector3(+b, +d, +c));
                    verts.Add(new Vector3(+b, +d, -c));
                    verts.Add(new Vector3(+b, -d, +c));
                    verts.Add(new Vector3(+b, -d, -c));

                    verts.Add(new Vector3(+d, +c, +b));
                    verts.Add(new Vector3(+d, -c, +b));
                    verts.Add(new Vector3(-d, +c, +b));
                    verts.Add(new Vector3(-d, -c, +b));

                    verts.Add(new Vector3(+c, +b, +d));
                    verts.Add(new Vector3(-c, +b, +d));
                    verts.Add(new Vector3(+c, +b, -d));
                    verts.Add(new Vector3(-c, +b, -d));
                    return verts.ToArray();
                }
            case PlatonicSolid.Square:
                return new Vector3[]
                {
                    new Vector3(-1, -1, 0),
                    new Vector3(-1,  1, 0),
                    new Vector3( 1,  1, 0),
                    new Vector3( 1, -1, 0)
                };
            case PlatonicSolid.Circle:
                // For Circle, determine the number of vertices from candidateSubdivision (clamped between 3 and spherePointCount).
                int circleVertexCount = Mathf.Clamp(candidateSubdivision, 3, spherePointCount);
                Vector3[] circleVerts = new Vector3[circleVertexCount];
                // Evenly distribute vertices around the circle.
                for (int i = 0; i < circleVertexCount; i++)
                {
                    float angle = 2 * Mathf.PI * i / circleVertexCount;
                    circleVerts[i] = new Vector3(Mathf.Cos(angle), Mathf.Sin(angle), 0);
                }
                return circleVerts;
            default:
                return null;
        }
    }

    /// <summary>
    /// Returns a list of faces for the selected solid.
    /// For 3D solids, each face is represented by an array of vertex indices.
    /// For 2D shapes, a single face (the entire polygon) is returned.
    /// </summary>
    /// <param name="solidType">The selected Platonic solid.</param>
    /// <returns>A list of integer arrays, each representing a face.</returns>
    private List<int[]> GetFaces(PlatonicSolid solidType)
    {
        List<int[]> faces = new List<int[]>();
        switch (solidType)
        {
            case PlatonicSolid.Tetrahedron:
                faces.Add(new int[] { 0, 1, 2 });
                faces.Add(new int[] { 0, 3, 1 });
                faces.Add(new int[] { 0, 2, 3 });
                faces.Add(new int[] { 1, 3, 2 });
                break;
            case PlatonicSolid.Cube:
                faces.Add(new int[] { 0, 1, 3, 2 });
                faces.Add(new int[] { 4, 6, 7, 5 });
                faces.Add(new int[] { 0, 4, 5, 1 });
                faces.Add(new int[] { 2, 3, 7, 6 });
                faces.Add(new int[] { 1, 5, 7, 3 });
                faces.Add(new int[] { 0, 2, 6, 4 });
                break;
            case PlatonicSolid.Octahedron:
                faces.Add(new int[] { 0, 2, 4 });
                faces.Add(new int[] { 2, 1, 4 });
                faces.Add(new int[] { 1, 3, 4 });
                faces.Add(new int[] { 3, 0, 4 });
                faces.Add(new int[] { 0, 5, 2 });
                faces.Add(new int[] { 2, 5, 1 });
                faces.Add(new int[] { 1, 5, 3 });
                faces.Add(new int[] { 3, 5, 0 });
                break;
            case PlatonicSolid.Icosahedron:
                faces.Add(new int[] { 0, 2, 8 });
                faces.Add(new int[] { 0, 4, 6 });
                faces.Add(new int[] { 0, 6, 10 });
                faces.Add(new int[] { 0, 8, 4 });
                faces.Add(new int[] { 0, 10, 2 });
                faces.Add(new int[] { 1, 3, 11 });
                faces.Add(new int[] { 1, 4, 9 });
                faces.Add(new int[] { 1, 6, 4 });
                faces.Add(new int[] { 1, 9, 3 });
                faces.Add(new int[] { 1, 11, 6 });
                faces.Add(new int[] { 2, 5, 8 });
                faces.Add(new int[] { 2, 7, 5 });
                faces.Add(new int[] { 2, 10, 7 });
                faces.Add(new int[] { 3, 5, 7 });
                faces.Add(new int[] { 3, 7, 11 });
                faces.Add(new int[] { 3, 9, 5 });
                faces.Add(new int[] { 4, 8, 9 });
                faces.Add(new int[] { 5, 9, 8 });
                faces.Add(new int[] { 6, 11, 10 });
                faces.Add(new int[] { 7, 10, 11 });
                break;
            case PlatonicSolid.Dodecahedron:
                faces.Add(new int[] { 0, 8, 10, 2, 16 });
                faces.Add(new int[] { 0, 16, 18, 1, 12 });
                faces.Add(new int[] { 0, 12, 14, 3, 8 });
                faces.Add(new int[] { 1, 9, 5, 14, 12 });
                faces.Add(new int[] { 1, 18, 4, 11, 9 });
                faces.Add(new int[] { 2, 10, 6, 15, 13 });
                faces.Add(new int[] { 2, 13, 4, 18, 16 });
                faces.Add(new int[] { 3, 14, 5, 19, 17 });
                faces.Add(new int[] { 3, 17, 6, 10, 8 });
                faces.Add(new int[] { 4, 13, 15, 7, 11 });
                faces.Add(new int[] { 5, 9, 11, 7, 19 });
                faces.Add(new int[] { 6, 17, 19, 7, 15 });
                break;
            case PlatonicSolid.Square:
                faces.Add(new int[] { 0, 1, 2, 3 });
                break;
            case PlatonicSolid.Circle:
                // For the Circle, determine the number of vertices based on candidateSubdivision (clamped).
                int circleVertexCount = Mathf.Clamp(candidateSubdivision, 3, spherePointCount);
                int[] indices = new int[circleVertexCount];
                for (int i = 0; i < circleVertexCount; i++)
                    indices[i] = i;
                faces.Add(indices);
                break;
        }
        return faces;
    }

    /// <summary>
    /// Makes a mesh double-sided by appending reversed triangles to the existing triangle list.
    /// This ensures that both the front and back faces are rendered.
    /// </summary>
    /// <param name="mesh">The mesh to modify.</param>
    /// <returns>The modified mesh with double-sided faces.</returns>
    private Mesh MakeMeshDoubleSided(Mesh mesh)
    {
        int[] originalTriangles = mesh.triangles;
        int triangleCount = originalTriangles.Length;
        int[] doubleTriangles = new int[triangleCount * 2];
        // Copy the original triangles.
        for (int i = 0; i < triangleCount; i++)
        {
            doubleTriangles[i] = originalTriangles[i];
        }
        // For each triangle, add a reversed triangle.
        for (int i = 0; i < triangleCount; i += 3)
        {
            doubleTriangles[triangleCount + i] = originalTriangles[i];
            doubleTriangles[triangleCount + i + 1] = originalTriangles[i + 2];
            doubleTriangles[triangleCount + i + 2] = originalTriangles[i + 1];
        }
        mesh.triangles = doubleTriangles;
        mesh.RecalculateNormals();
        return mesh;
    }


    #endregion
}

#endregion
