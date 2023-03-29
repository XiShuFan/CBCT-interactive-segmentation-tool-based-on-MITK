project(org_mitk_panorama)

mitk_create_plugin(
    EXPORT_DIRECTIVE Panorama_EXPORT
	EXPORTED_INCLUDE_SUFFIXES src
	MODULE_DEPENDS MitkQtWidgetsExt MitkTestingHelper MitkModelFit MitkModelFitUI
)
