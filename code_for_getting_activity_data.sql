-- SELECT assay_data->'sid' AS sid FROM activity;
SELECT json_array_elements(assay_data) FROM activity WHERE bioassay_id = '382001';