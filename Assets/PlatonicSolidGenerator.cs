#if UNITY_EDITOR
using UnityEditor;
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
    Tetrahedron,
    Cube,
    Octahedron,
    Icosahedron,
    Dodecahedron,
    Square,   // 2D: a single face (from a cube)
    Circle    // 2D: a circle with evenly spaced vertices
}

/// <summary>
/// Combined connection ordering enum. The first five values use legacy coordinate‑ordering;
/// the remaining four produce a spiral ordering.
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
/// PlatonicSolidGizmo draws a wireframe of a chosen Platonic (or pseudo‑Platonic) solid,
/// generates a cloud of sphere points over its faces, assigns letters to each sphere,
/// processes an active letter chain (which can be obfuscated via several string manipulation methods)
/// to filter and connect a chain of spheres, and draws connecting cylinders between them.
/// When the "drawEndpointsOnly" option is enabled, only the first and last spheres are drawn
/// (with the last replaced by a cube), but labels are drawn at every active point.
/// A label below the solid (if enabled) displays the processed active filter.
/// </summary>
public class PlatonicSolidGizmo : MonoBehaviour
{
    #region Public Variables

    [Header("Solid Settings")]
    [Tooltip("Choose which Platonic solid to display")]
    public PlatonicSolid solid = PlatonicSolid.Tetrahedron;
    [Tooltip("Uniform scale for the solid")]
    public float scale = 1f;

    [Header("Wireframe Options")]
    [Tooltip("Color of the primary wireframe of the Platonic solid")]
    public Color wireframeColor = Color.green;
    [Tooltip("Toggle primary wireframe drawing on or off")]
    public bool drawWireframe = true;

    [Header("Sphere Cloud Settings")]
    [Tooltip("Exact total number of spheres to distribute over the faces")]
    public int spherePointCount = 50;
    [Tooltip("Base radius for the first sphere")]
    public float sphereRadius = 0.05f;
    [Tooltip("Additional radius added per sphere (each subsequent sphere is larger)")]
    public float sphereSizeIncrement = 0.01f;

    [Header("Sphere Color Gradient")]
    [Tooltip("Color (with opacity) for the first sphere")]
    public Color startSphereColor = new Color(1f, 0f, 0f, 1f);
    [Tooltip("Color (with opacity) for the last sphere")]
    public Color endSphereColor = new Color(0f, 0f, 1f, 0.5f);

    [Header("Letter Mapping Settings")]
    [Tooltip("If false, letters are assigned in A–Z order; if true, use alternate ordering")]
    public bool useAlternateAlphabet = false;
    [Tooltip("Alphabet used for A–Z mapping")]
    public string alphabetAZ = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    [Tooltip("Alternate alphabet ordering")]
    public string alphabetEQ = "ESIARNTOLCDUGPMHBYFVKWZXJQ";
    [Tooltip("Active letter chain. Type a sentence here which will be processed (via obfuscation methods) " +
             "to produce the chain of active letters. For example, 'AGPA' selects spheres with A, then G, then P, then A (even if letters repeat)")]
    public string activeLettersFilter = "";

    [Header("Active Letters Obfuscation Settings")]
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

    [Header("Endpoint Options")]
    [Tooltip("When enabled, only the first and last spheres are drawn (the last is drawn as a cube). " +
             "Intermediate sphere points are not drawn, but labels are still placed at these points.")]
    public bool drawEndpointsOnly = false;

    [Header("Label Options")]
    [Tooltip("Toggle all labels (sphere labels, processed filter label, etc.) on or off")]
    public bool drawLabels = true;

    [Header("Candidate Generation Settings")]
    [Tooltip("Subdivision level for candidate point generation on each face (for Circle, clamped to determine vertex count)")]
    public int candidateSubdivision = 10;

    [Header("Cylinder Settings")]
    [Tooltip("Radius for the connecting cylinders between sphere centers")]
    public float cylinderRadius = 0.05f;

    [Header("Connection Ordering Settings")]
    [Tooltip("Combined connection setting that determines how sphere points are ordered. " +
             "Legacy options (TopToBottom, etc.) sort by coordinate; spiral options (ClockwiseOut/In, CounterclockwiseOut/In) produce a spiral." +
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
    /// Called by Unity to render the gizmos. This method draws (if enabled) the primary wireframe,
    /// generates sphere points over the solid, orders them (using legacy or spiral ordering),
    /// maps letters to each sphere, processes and applies the active letter chain (using obfuscation methods),
    /// and draws the visible spheres (with labels if enabled) and connecting cylinders.
    /// When drawEndpointsOnly is enabled, only the first and last spheres are drawn (with the last replaced by a cube),
    /// but labels are drawn at every active point.
    /// A label below the solid (if enabled) displays the processed active filter.
    /// </summary>
    private void OnDrawGizmos()
    {
        // 1. Retrieve and transform vertices for the primary solid (scaled by "scale").
        Vector3[] vertices = GetVertices(solid);
        if (vertices == null || vertices.Length == 0)
            return;
        for (int i = 0; i < vertices.Length; i++)
            vertices[i] = transform.position + vertices[i] * scale;

        // 2. Optionally draw the primary wireframe.
        if (drawWireframe)
        {
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

        // 3. Generate the final set of sphere points.
        List<Vector3> finalPoints = new List<Vector3>();
        if (spherePointCount == vertices.Length)
        {
            finalPoints.AddRange(vertices);
        }
        else
        {
            List<int[]> faces = GetFaces(solid);
            if (faces == null || faces.Count == 0)
                return;
            List<Vector3> candidatePoints = new List<Vector3>();
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

        // 4. Order the sphere points according to the connection setting.
        List<Vector3> orderedPoints = new List<Vector3>();
        if (connectionSetting == ConnectionSetting.TopToBottom ||
            connectionSetting == ConnectionSetting.BottomToTop ||
            connectionSetting == ConnectionSetting.LeftToRight ||
            connectionSetting == ConnectionSetting.RightToLeft ||
            connectionSetting == ConnectionSetting.Random)
        {
            // Legacy ordering: coordinate‐based.
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
                    for (int i = 0; i < orderedPoints.Count; i++)
                    {
                        int rnd = UnityEngine.Random.Range(i, orderedPoints.Count);
                        Vector3 temp = orderedPoints[i];
                        orderedPoints[i] = orderedPoints[rnd];
                        orderedPoints[rnd] = temp;
                    }
                    break;
            }
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
            // Spiral ordering: use polar coordinates combined with a weighted proximity term.
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

        // 5. Map letters to each ordered point.
        int pointCount = orderedPoints.Count;
        string[] letterMapping = new string[pointCount];
        for (int i = 0; i < pointCount; i++)
            letterMapping[i] = "";
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

        // 6. Process the activeLettersFilter using obfuscation methods.
        string processedFilter = ProcessActiveLettersFilter(activeLettersFilter);

        // 7. Build the active chain from the processed filter.
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

        // 8. Build lists of visible points and assign gradient colors.
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

        // 9. Draw the visible spheres (or endpoints) and draw labels at every active point.
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

        // 10. Draw connecting cylinders between consecutive visible points.
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

        // 11. Draw a centered label below the solid displaying the processed active filter.
#if UNITY_EDITOR
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
    /// The methods are applied in the following order:
    /// ToUpper (always), then optionally: Flip, Rotate, CutUp, RemoveWhitespace, Reverse, Shift, Disemvowel, Squeeze, Deduplicate, and Rearrange.
    /// </summary>
    /// <param name="filter">The original filter string.</param>
    /// <returns>The processed (obfuscated) string.</returns>
    private string ProcessActiveLettersFilter(string filter)
    {
        if (string.IsNullOrEmpty(filter))
            return filter;
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
    /// Standard farthest point sampling. Selects exactly 'count' points from the candidate list,
    /// starting with a candidate that is near one of the original vertices.
    /// </summary>
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
    /// Modified farthest point sampling that always includes a preselected set of points (e.g. vertices)
    /// and fills the remaining points from the candidate set.
    /// </summary>
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
    /// Generates candidate points on a polygon by fan‑triangulating the polygon and creating a barycentric grid.
    /// </summary>
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
    /// Removes duplicate points from the list using the specified tolerance.
    /// </summary>
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
    /// Returns the vertices for the chosen Platonic (or pseudo‑Platonic) solid.
    /// 3D solids are defined using symmetric coordinates; Square and Circle are defined in the XY plane.
    /// </summary>
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
                float phi = (1f + Mathf.Sqrt(5f)) / 2f;
                return new Vector3[]
                {
                    new Vector3(0,  1,  phi),
                    new Vector3(0,  1, -phi),
                    new Vector3(0, -1,  phi),
                    new Vector3(0, -1, -phi),
                    new Vector3( 1,  phi, 0),
                    new Vector3( 1, -phi, 0),
                    new Vector3(-1,  phi, 0),
                    new Vector3(-1, -phi, 0),
                    new Vector3( phi, 0,  1),
                    new Vector3( phi, 0, -1),
                    new Vector3(-phi, 0,  1),
                    new Vector3(-phi, 0, -1)
                };
            case PlatonicSolid.Dodecahedron:
                float phiD = (1f + Mathf.Sqrt(5f)) / 2f;
                float invPhi = 1f / phiD;
                List<Vector3> verts = new List<Vector3>();
                verts.Add(new Vector3(1, 1, 1));
                verts.Add(new Vector3(1, 1, -1));
                verts.Add(new Vector3(1, -1, 1));
                verts.Add(new Vector3(1, -1, -1));
                verts.Add(new Vector3(-1, 1, 1));
                verts.Add(new Vector3(-1, 1, -1));
                verts.Add(new Vector3(-1, -1, 1));
                verts.Add(new Vector3(-1, -1, -1));
                verts.Add(new Vector3(0, invPhi, phiD));
                verts.Add(new Vector3(0, invPhi, -phiD));
                verts.Add(new Vector3(0, -invPhi, phiD));
                verts.Add(new Vector3(0, -invPhi, -phiD));
                verts.Add(new Vector3(invPhi, phiD, 0));
                verts.Add(new Vector3(invPhi, -phiD, 0));
                verts.Add(new Vector3(-invPhi, phiD, 0));
                verts.Add(new Vector3(-invPhi, -phiD, 0));
                verts.Add(new Vector3(phiD, 0, invPhi));
                verts.Add(new Vector3(phiD, 0, -invPhi));
                verts.Add(new Vector3(-phiD, 0, invPhi));
                verts.Add(new Vector3(-phiD, 0, -invPhi));
                return verts.ToArray();
            case PlatonicSolid.Square:
                return new Vector3[]
                {
                    new Vector3(-1, -1, 0),
                    new Vector3(-1,  1, 0),
                    new Vector3( 1,  1, 0),
                    new Vector3( 1, -1, 0)
                };
            case PlatonicSolid.Circle:
                int circleVertexCount = Mathf.Clamp(candidateSubdivision, 3, spherePointCount);
                Vector3[] circleVerts = new Vector3[circleVertexCount];
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
    /// Returns a list of faces for the chosen solid.
    /// For 3D solids, each face is represented as an array of vertex indices.
    /// For 2D shapes, a single face (the polygon) is returned.
    /// </summary>
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
                faces.Add(new int[] { 0, 8, 4 });
                faces.Add(new int[] { 0, 4, 6 });
                faces.Add(new int[] { 0, 6, 10 });
                faces.Add(new int[] { 0, 10, 2 });
                faces.Add(new int[] { 0, 2, 8 });
                faces.Add(new int[] { 8, 2, 5 });
                faces.Add(new int[] { 8, 5, 9 });
                faces.Add(new int[] { 8, 9, 4 });
                faces.Add(new int[] { 4, 9, 1 });
                faces.Add(new int[] { 4, 1, 6 });
                faces.Add(new int[] { 6, 1, 11 });
                faces.Add(new int[] { 6, 11, 10 });
                faces.Add(new int[] { 10, 11, 3 });
                faces.Add(new int[] { 10, 3, 2 });
                faces.Add(new int[] { 2, 3, 5 });
                faces.Add(new int[] { 5, 3, 7 });
                faces.Add(new int[] { 5, 7, 9 });
                faces.Add(new int[] { 9, 7, 1 });
                faces.Add(new int[] { 1, 7, 11 });
                faces.Add(new int[] { 11, 7, 3 });
                break;
            case PlatonicSolid.Dodecahedron:
                faces.Add(new int[] { 0, 16, 2, 10, 8 });
                faces.Add(new int[] { 0, 8, 4, 14, 12 });
                faces.Add(new int[] { 0, 12, 1, 17, 16 });
                faces.Add(new int[] { 1, 12, 14, 5, 19 });
                faces.Add(new int[] { 1, 19, 7, 11, 17 });
                faces.Add(new int[] { 2, 16, 17, 11, 6 });
                faces.Add(new int[] { 2, 10, 3, 13, 6 });
                faces.Add(new int[] { 3, 10, 8, 9, 13 });
                faces.Add(new int[] { 4, 8, 9, 15, 14 });
                faces.Add(new int[] { 5, 14, 15, 7, 19 });
                faces.Add(new int[] { 6, 13, 9, 11, 7 });
                faces.Add(new int[] { 17, 11, 9, 15, 19 });
                break;
            case PlatonicSolid.Square:
                faces.Add(new int[] { 0, 1, 2, 3 });
                break;
            case PlatonicSolid.Circle:
                int circleVertexCount = Mathf.Clamp(candidateSubdivision, 3, spherePointCount);
                int[] indices = new int[circleVertexCount];
                for (int i = 0; i < circleVertexCount; i++)
                    indices[i] = i;
                faces.Add(indices);
                break;
        }
        return faces;
    }

    #endregion
}

#endregion

#region String Obfuscation Extensions

/// <summary>
/// A collection of extension methods to obfuscate a string. These methods are applied
/// (in a fixed order) to the activeLettersFilter string so that you can type a sentence
/// and have it obfuscated before mapping its letters to sphere points.
/// </summary>
public static class StringExtensions
{
    public static string Disemvowel(this string str)
    {
        return new string(str.Where(c => !"aeiouAEIOU".Contains(c)).ToArray());
    }

    public static string Squeeze(this string str)
    {
        StringBuilder sb = new StringBuilder();
        foreach (char c in str)
        {
            if (sb.Length == 0 || sb[sb.Length - 1] != c)
                sb.Append(c);
        }
        return sb.ToString();
    }

    public static string Deduplicate(this string str)
    {
        return new string(str.Distinct().ToArray());
    }

    public static string RemoveWhitespace(this string str)
    {
        StringBuilder sb = new StringBuilder();
        foreach (char c in str)
        {
            if (!char.IsWhiteSpace(c))
                sb.Append(c);
        }
        return sb.ToString();
    }

    public static string Rearrange(this string str)
    {
        char[] chars = str.ToCharArray();
        System.Random rng = new System.Random();
        for (int i = chars.Length - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (chars[i], chars[j]) = (chars[j], chars[i]);
        }
        return new string(chars);
    }

    public static string Shift(this string str, int? offset = null)
    {
        int actualOffset = offset ?? new System.Random().Next(26);
        StringBuilder sb = new StringBuilder();
        foreach (char c in str)
        {
            if (char.IsLetter(c))
            {
                char baseChar = char.IsUpper(c) ? 'A' : 'a';
                int shifted = (c - baseChar + actualOffset) % 26;
                sb.Append((char)(baseChar + shifted));
            }
            else
                sb.Append(c);
        }
        return sb.ToString();
    }

    public static string Reverse(this string str)
    {
        char[] arr = str.ToCharArray();
        Array.Reverse(arr);
        return new string(arr);
    }

    public static string CutUp(this string str)
    {
        string[] words = str.Split(' ');
        System.Random rng = new System.Random();
        for (int i = words.Length - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (words[i], words[j]) = (words[j], words[i]);
        }
        return string.Join(" ", words);
    }

    public static string Rotate(this string str, int? offset = null)
    {
        string[] words = str.Split(' ');
        int actualOffset = offset ?? new System.Random().Next(words.Length);
        string[] rotated = new string[words.Length];
        for (int i = 0; i < words.Length; i++)
            rotated[(i + actualOffset) % words.Length] = words[i];
        return string.Join(" ", rotated);
    }

    public static string Flip(this string str)
    {
        string[] words = str.Split(' ');
        Array.Reverse(words);
        return string.Join(" ", words);
    }
}

#endregion
